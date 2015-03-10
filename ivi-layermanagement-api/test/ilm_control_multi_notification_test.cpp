/***************************************************************************
 *
 * Copyright 2015 Codethink Ltd
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 ****************************************************************************/

#include <gtest/gtest.h>
#include <stdio.h>

#include <unistd.h>
#include <sys/types.h>
#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>
#include <errno.h>
#include <stdlib.h>
#include <signal.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>


#include "TestBase.h"

extern "C" {
    #include "ilm_client.h"
    #include "ilm_control.h"
}

void add_n_secs(struct timespec *tv, long nsec)
{
   assert(nsec < 1000000000);

   tv->tv_nsec += nsec;
   if (tv->tv_nsec >= 1000000000)
   {
      tv->tv_nsec -= 1000000000;
      tv->tv_sec++;
   }
}

template <typename T>
bool contains(T const *actual, size_t as, T expected)
{
   for (unsigned i = 0; i < as; i++)
      if (actual[i] == expected)
         return true;
   return false;
}

/* Tests with callbacks
 * For each test first set the global variables to point to where parameters of the callbacks are supposed to be placed.
 */
static pthread_mutex_t notificationMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  waiterVariable = PTHREAD_COND_INITIALIZER;
static int timesCalled=0;

struct surface_def {
    t_ilm_uint requestedSurfaceId;
    t_ilm_uint returnedSurfaceId;
    ilmSurfaceProperties surfaceProperties;
};

struct layer_def {
    t_ilm_uint layerId;
    ilmLayerProperties layerProperties;
};

class IlmMultiNotificationTest : public TestBase, public ::testing::Test {

    class PthreadMutexLock {
    private:
        pthread_mutex_t &mutex_;

        PthreadMutexLock(const PthreadMutexLock&);
        PthreadMutexLock& operator=(const PthreadMutexLock&);

    public:
        explicit PthreadMutexLock(pthread_mutex_t &mutex): mutex_(mutex) {
            pthread_mutex_lock(&mutex_);
        }
        ~PthreadMutexLock() {
            pthread_mutex_unlock(&mutex_);
        }
    };

public:

    t_ilm_uint layer;

    // Pointers where to put received values for current Test
    static t_ilm_layer callbackLayerId;
    static t_ilm_surface callbackSurfaceId;
    static struct ilmLayerProperties LayerProperties;
    static unsigned int mask;
    static t_ilm_surface surface;
    static ilmSurfaceProperties SurfaceProperties;

    std::vector<layer_def> layers_allocated;
    std::vector<surface_def> surfaces_allocated;

    IlmMultiNotificationTest(){}
    ~IlmMultiNotificationTest(){}

    void SetUp()
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_initWithNativedisplay((t_ilm_nativedisplay)wlDisplay));
        mask = static_cast<t_ilm_notification_mask>(0);
        surface = -1;
        timesCalled=0;

        callbackLayerId = INVALID_ID;
        callbackSurfaceId = INVALID_ID;
    }

    void TearDown()
    {
        // set default values
        callbackLayerId = -1;
        t_ilm_layer* layers = NULL;
        t_ilm_int numLayer=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&numLayer, &layers));
        for (t_ilm_int i=0; i<numLayer; i++)
        {
            if ( ( layers[i] >= getStartLayerId() ) && ( layers[i] <= getEndLayerId() ) )
            {
                EXPECT_EQ( ILM_SUCCESS, ilm_layerRemove( layers[i] ) );
            }
        };

        layers_allocated.clear();

        free(layers);

        t_ilm_surface* surfaces = NULL;
        t_ilm_int numSurfaces=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&numSurfaces, &surfaces));
        for (t_ilm_int i=0; i<numSurfaces; i++)
        {
            if ( ( surfaces[i] >= getStartSurfaceId() ) && ( layers[i] <= getEndSurfaceId() ) )
            {
                EXPECT_EQ( ILM_SUCCESS, ilm_surfaceRemove( surfaces[i] ) );
            }
        };

        free(surfaces);

        surfaces_allocated.clear();

        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());
        EXPECT_EQ(ILM_SUCCESS, ilm_destroy());
    }


    static void assertCallbackcalled(int numberOfExpectedCalls=1){
        static struct timespec theTime;
        clock_gettime(CLOCK_REALTIME, &theTime);
        add_n_secs(&theTime, 500000000);
        PthreadMutexLock lock(notificationMutex);
        int status = 0;
        do {
            if (numberOfExpectedCalls!=timesCalled) {
                status = pthread_cond_timedwait( &waiterVariable, &notificationMutex, &theTime);
            }
        } while (status!=ETIMEDOUT && numberOfExpectedCalls!=timesCalled);

        // we cannot rely on a timeout as layer callbacks are always called synchronously on ilm_commitChanges()
        EXPECT_NE(ETIMEDOUT, status);
        timesCalled=0;
    }

    static void assertNoCallbackIsCalled(){
        struct timespec theTime;
        clock_gettime(CLOCK_REALTIME, &theTime);
        add_n_secs(&theTime, 500000000);
        PthreadMutexLock lock(notificationMutex);
        // assert that we have not been notified
        ASSERT_EQ(ETIMEDOUT, pthread_cond_timedwait( &waiterVariable, &notificationMutex, &theTime));
    }

    static void LayerCallbackFunction(t_ilm_layer layer, struct ilmLayerProperties* layerProperties, t_ilm_notification_mask m)
    {
        PthreadMutexLock lock(notificationMutex);

        if ((unsigned)m & ILM_NOTIFICATION_VISIBILITY)
        {
            LayerProperties.visibility = layerProperties->visibility;
        }

        if ((unsigned)m & ILM_NOTIFICATION_OPACITY)
        {
            LayerProperties.opacity = layerProperties->opacity;
        }

        if ((unsigned)m & ILM_NOTIFICATION_ORIENTATION)
        {
            LayerProperties.orientation = layerProperties->orientation;
        }

        if ((unsigned)m & ILM_NOTIFICATION_SOURCE_RECT)
        {
            LayerProperties.sourceX = layerProperties->sourceX;
            LayerProperties.sourceY = layerProperties->sourceY;
            LayerProperties.sourceWidth = layerProperties->sourceWidth;
            LayerProperties.sourceHeight = layerProperties->sourceHeight;
        }

        if ((unsigned)m & ILM_NOTIFICATION_DEST_RECT)
        {
            LayerProperties.destX = layerProperties->destX;
            LayerProperties.destY = layerProperties->destY;
            LayerProperties.destWidth = layerProperties->destWidth;
            LayerProperties.destHeight = layerProperties->destHeight;
        }

        EXPECT_TRUE(callbackLayerId == (unsigned)-1 || callbackLayerId == layer);
        callbackLayerId = layer;
        mask |= (unsigned)m;
        timesCalled++;

        pthread_cond_signal( &waiterVariable );
    }

    static void SurfaceCallbackFunction(t_ilm_surface surface, struct ilmSurfaceProperties* surfaceProperties, t_ilm_notification_mask m)
    {
        PthreadMutexLock lock(notificationMutex);

        if ((unsigned)m & ILM_NOTIFICATION_VISIBILITY)
        {
            SurfaceProperties.visibility = surfaceProperties->visibility;
        }

        if ((unsigned)m & ILM_NOTIFICATION_OPACITY)
        {
            SurfaceProperties.opacity = surfaceProperties->opacity;
        }

        if ((unsigned)m & ILM_NOTIFICATION_ORIENTATION)
        {
            SurfaceProperties.orientation = surfaceProperties->orientation;
        }

        if ((unsigned)m & ILM_NOTIFICATION_SOURCE_RECT)
        {
            SurfaceProperties.sourceX = surfaceProperties->sourceX;
            SurfaceProperties.sourceY = surfaceProperties->sourceY;
            SurfaceProperties.sourceWidth = surfaceProperties->sourceWidth;
            SurfaceProperties.sourceHeight = surfaceProperties->sourceHeight;
        }

        if ((unsigned)m & ILM_NOTIFICATION_DEST_RECT)
        {
            SurfaceProperties.destX = surfaceProperties->destX;
            SurfaceProperties.destY = surfaceProperties->destY;
            SurfaceProperties.destWidth = surfaceProperties->destWidth;
            SurfaceProperties.destHeight = surfaceProperties->destHeight;
        }

//        EXPECT_TRUE(callbackSurfaceId == (unsigned)-1 || callbackSurfaceId == surface);
        callbackSurfaceId = surface;
        mask |= (unsigned)m;
        timesCalled++;

        pthread_cond_signal( &waiterVariable );
    }
};

// Pointers where to put received values for current Test
t_ilm_layer IlmMultiNotificationTest::callbackLayerId;
t_ilm_surface IlmMultiNotificationTest::callbackSurfaceId;
struct ilmLayerProperties IlmMultiNotificationTest::LayerProperties;
unsigned int IlmMultiNotificationTest::mask;
t_ilm_surface IlmMultiNotificationTest::surface;
ilmSurfaceProperties IlmMultiNotificationTest::SurfaceProperties;

TEST_F(IlmMultiNotificationTest, ilm_multiLayerAddNotificationWithoutCallback)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 15 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 25 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = surf_pos[0];
        layers_allocated[i].layerProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect no callback to have been called
        assertNoCallbackIsCalled();
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect no callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {

            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (IDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                // Iterate round remaining layers and check dimensions 
                for (int k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                    // Change something that has been pre-set and check callback
                    ASSERT_EQ(ILM_SUCCESS, ilm_layerSetVisibility(layers_allocated[index].layerId, ILM_TRUE));

                    // expect callback to have been called
                    assertNoCallbackIsCalled();

                }
            }
        }

        free(IDs);
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (IDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }

        free(IDs);
    }

    surfaces_allocated.clear();


}

TEST_F(IlmMultiNotificationTest, ilm_multiSurfaceAddNotificationWithoutCallback)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    e_ilmPixelFormat pixelFormats[8] = {
        ILM_PIXELFORMAT_R_8,
        ILM_PIXELFORMAT_RGB_888,
        ILM_PIXELFORMAT_RGBA_8888,
        ILM_PIXELFORMAT_RGB_565,
        ILM_PIXELFORMAT_RGBA_5551,
        ILM_PIXELFORMAT_RGBA_6661,
        ILM_PIXELFORMAT_RGBA_4444,
        ILM_PIXEL_FORMAT_UNKNOWN };


    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 15 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 25 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     pixelFormats[i],
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Add notifications to surface
    for (uint i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceAddNotification(surfaces_allocated[i].returnedSurfaceId,
                  &SurfaceCallbackFunction));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces check call backs
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_pos[2] = {22 + (i * 5), 42 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        callbackSurfaceId = surfaces_allocated[i].returnedSurfaceId;                
        assertCallbackcalled();
        EXPECT_EQ(callbackSurfaceId,
                  surfaces_allocated[i].returnedSurfaceId);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm parameters are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (IDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceX,
                          posreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceY,
                          posreturned[1]);

                // Change position of surfaces and check callback
                t_ilm_uint pos_set[2] = {23 + (i * 5), 44 + (i * 5)};
                surfaces_allocated[index].surfaceProperties.sourceX = pos_set[0];
                surfaces_allocated[index].surfaceProperties.sourceY = pos_set[1];
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_surfaceSetPosition(surfaces_allocated[index].returnedSurfaceId,
                          pos_set));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
                callbackSurfaceId = surfaces_allocated[index].returnedSurfaceId;                
                assertCallbackcalled();
                EXPECT_EQ(callbackSurfaceId,
                          surfaces_allocated[index].returnedSurfaceId);

            }
        }

        free(IDs);
    }

    surfaces_allocated.clear();
}

TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetPosition)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 15 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 25 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // Set position of layers check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        callbackLayerId = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have not been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemoveNotification(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                // Iterate round remaining layers and check dimensions 
                for (int k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    callbackLayerId = layerIDs[j];

                    EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                    // Change something that has been pre-set and check callback
                    ASSERT_EQ(ILM_SUCCESS, ilm_layerSetVisibility(layerIDs[j], ILM_TRUE));
                    ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                    assertCallbackcalled();
                    EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);

                }
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }
    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetDimension)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 15 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 25 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Change dimensions of layers check callbacks
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_dim[2] = {30 + (i * 5), 20 + (i * 5)};
        callbackLayerId = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.origSourceWidth = layer_dim[0];
        layers_allocated[i].layerProperties.origSourceHeight = layer_dim[1];

        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetDimension(layers_allocated[i].layerId, layer_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);
    }

    // Set position of layers check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        callbackLayerId = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.sourceX = surf_pos[0];
        layers_allocated[i].layerProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);

    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to not have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                // Iterate round remaining layers and check dimensions 
                for (int k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    callbackLayerId = layers_allocated[index].layerId;

                    EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                    // Change something that has been pre-set and check callback
                    ASSERT_EQ(ILM_SUCCESS, ilm_layerSetVisibility(layers_allocated[index].layerId, ILM_TRUE));
                    EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

                    // expect callback to have been called
                    assertCallbackcalled();
                    EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);

                }
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}

TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetDestinationRectangle)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 12 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 32 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {30 + (i * 5), 20 + (i * 5)};
        callbackLayerId = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.destX = 3 + (i * 5);
        layers_allocated[i].layerProperties.destY = 2 + (i * 5);
        layers_allocated[i].layerProperties.destWidth = 50 + (i * 10);
        layers_allocated[i].layerProperties.destHeight = 55 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetDestinationRectangle(layers_allocated[i].layerId,
                                                   layers_allocated[i].layerProperties.destX,
                                                   layers_allocated[i].layerProperties.destY,
                                                   layers_allocated[i].layerProperties.destWidth,
                                                   layers_allocated[i].layerProperties.destHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.destWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.destHeight, dimreturned[1]);

                // Confirm destination rectangle
                ilmLayerProperties layerProperties;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                ASSERT_EQ(layers_allocated[index].layerProperties.destX,
                          layerProperties.destX);
                ASSERT_EQ(layers_allocated[index].layerProperties.destY,
                          layerProperties.destY);
                ASSERT_EQ(layers_allocated[index].layerProperties.destWidth,
                          layerProperties.destWidth);
                ASSERT_EQ(layers_allocated[index].layerProperties.destHeight,
                          layerProperties.destHeight);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetOpacity)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 18 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 39 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 202 * (i + 1);
        layer->layerProperties.origSourceHeight = 244 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set Opacity of surfaces
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.opacity = 0.05 + (i * 0.15);
        callbackLayerId = layers_allocated[i].layerId;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOpacity(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                // Confirm opacity
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layerIDs[j], &returned));
                EXPECT_NEAR(layers_allocated[index].layerProperties.opacity,
                            returned, 0.01);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetOrientation)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 22 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 42 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 202 * (i + 1);
        layer->layerProperties.origSourceHeight = 244 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = surf_pos[0];
        layers_allocated[i].layerProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        layers_allocated[i].layerProperties.orientation = orientations[no_layers % 4];
        callbackLayerId = layers_allocated[i].layerId;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetOrientation(layerIDs[j], &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetSourceRectangle)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 14 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 31 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {30 + (i * 5), 20 + (i * 5)};
        callbackLayerId = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.sourceX = 3 + (i * 5);
        layers_allocated[i].layerProperties.sourceY = 2 + (i * 5);
        layers_allocated[i].layerProperties.sourceWidth = 50 + (i * 10);
        layers_allocated[i].layerProperties.sourceHeight = 55 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;

                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                // Confirm source rectangle
                ilmLayerProperties layerProperties;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceX,
                          layerProperties.sourceX);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceY,
                          layerProperties.sourceY);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceWidth,
                          layerProperties.sourceWidth);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceHeight,
                          layerProperties.sourceHeight);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerSetVisibility)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 22 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 42 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 202 * (i + 1);
        layer->layerProperties.origSourceHeight = 244 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = surf_pos[0];
        layers_allocated[i].layerProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set Visibility
    for (int i = 0; i < no_layers; i++)
    {
        layers_allocated[i].layerProperties.visibility = visibility[no_layers % 2];
        callbackLayerId = layers_allocated[i].layerId;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                // Confirm visibility
                t_ilm_bool visibility_rtn;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetVisibility(layerIDs[j], &visibility_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.visibility,
                            visibility_rtn);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetOrientation(layers_allocated[index].layerId,
                          ILM_ZERO));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerMultipleValues1)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 16 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 28 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 198 * (i + 1);
        layer->layerProperties.origSourceHeight = 232 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = { 22 + (i * 5), 42 + (i * 5) };
        layers_allocated[i].layerProperties.sourceX = surf_pos[0];
        layers_allocated[i].layerProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;
        callbackLayerId = layers_allocated[i].layerId;
        t_ilm_uint surf_dim[2] = { 30 + (i * 5), 20 + (i * 5) };

        layers_allocated[i].layerProperties.sourceX = 3 + (i * 5);
        layers_allocated[i].layerProperties.sourceY = 2 + (i * 5);
        layers_allocated[i].layerProperties.sourceWidth = 50 + (i * 10);
        layers_allocated[i].layerProperties.sourceHeight = 55 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);
    }

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        callbackLayerId = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.orientation = orientations[i % 4];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                // Confirm destination rectangle
                ilmLayerProperties layerProperties;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceX,
                          layerProperties.sourceX);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceY,
                          layerProperties.sourceY);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceWidth,
                          layerProperties.sourceWidth);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceHeight,
                          layerProperties.sourceHeight);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_layerGetOrientation(layers_allocated[index].layerId,
                          &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Change something that has been pre-set and check callback
                layer = layers_allocated[index].layerId;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}

TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerMultipleValues2)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 26 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 48 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 198 * (i + 1);
        layer->layerProperties.origSourceHeight = 232 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = { 22 + (i * 5), 42 + (i * 5) };
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        callbackLayerId = layers_allocated[i].layerId;
        t_ilm_uint surf_dim[2] = { 30 + (i * 5), 20 + (i * 5) };

        layers_allocated[i].layerProperties.sourceX = 3 + (i * 5);
        layers_allocated[i].layerProperties.sourceY = 2 + (i * 5);
        layers_allocated[i].layerProperties.sourceWidth = 50 + (i * 10);
        layers_allocated[i].layerProperties.sourceHeight = 55 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId,callbackLayerId);
    }

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        callbackLayerId = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.orientation = orientations[i % 4];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Set Opacity of layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.opacity = 0.23 + (i * 0.15);
        callbackLayerId = layers_allocated[i].layerId;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOpacity(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            layer = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);

                // Confirm destination rectangle
                ilmLayerProperties layerProperties;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceX,
                          layerProperties.sourceX);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceY,
                          layerProperties.sourceY);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceWidth,
                          layerProperties.sourceWidth);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceHeight,
                          layerProperties.sourceHeight);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_layerGetOrientation(layers_allocated[index].layerId,
                          &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Confirm opacity
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layerIDs[j], &returned));
                EXPECT_NEAR(layers_allocated[index].layerProperties.opacity,
                            returned, 0.01);

                // Change something that has been pre-set and check callback
                layer = layers_allocated[index].layerId;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnLayerAllValues)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 28 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 49 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 198 * (i + 1);
        layer->layerProperties.origSourceHeight = 232 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = { 22 + (i * 5), 42 + (i * 5) };
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set source rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.sourceX = 3 + (i * 5);
        layers_allocated[i].layerProperties.sourceY = 2 + (i * 5);
        layers_allocated[i].layerProperties.sourceWidth = 50 + (i * 10);
        layers_allocated[i].layerProperties.sourceHeight = 55 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.destX = 5 + (i * 5);
        layers_allocated[i].layerProperties.destY = 3 + (i * 5);
        layers_allocated[i].layerProperties.destWidth = 51 + (i * 10);
        layers_allocated[i].layerProperties.destHeight = 63 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetDestinationRectangle(layers_allocated[i].layerId,
                                                   layers_allocated[i].layerProperties.destX,
                                                   layers_allocated[i].layerProperties.destY,
                                                   layers_allocated[i].layerProperties.destWidth,
                                                   layers_allocated[i].layerProperties.destHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        layer = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.orientation = orientations[i % 4];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set Opacity of layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.opacity = 0.23 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOpacity(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set Visibility
    for (int i = 0; i < no_layers; i++)
    {
        layers_allocated[i].layerProperties.visibility = visibility[no_layers % 2];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            layer = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.destWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.destHeight, dimreturned[1]);

                // Confirm source rectangle
                ilmLayerProperties layerProperties;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceX,
                          layerProperties.sourceX);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceY,
                          layerProperties.sourceY);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceWidth,
                          layerProperties.sourceWidth);
                ASSERT_EQ(layers_allocated[index].layerProperties.sourceHeight,
                          layerProperties.sourceHeight);

                // Confirm destination rectangle
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                ASSERT_EQ(layers_allocated[index].layerProperties.destX,
                          layerProperties.destX);
                ASSERT_EQ(layers_allocated[index].layerProperties.destY,
                          layerProperties.destY);
                ASSERT_EQ(layers_allocated[index].layerProperties.destWidth,
                          layerProperties.destWidth);
                ASSERT_EQ(layers_allocated[index].layerProperties.destHeight,
                          layerProperties.destHeight);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_layerGetOrientation(layers_allocated[index].layerId,
                          &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Confirm opacity
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layerIDs[j], &returned));
                EXPECT_NEAR(layers_allocated[index].layerProperties.opacity,
                            returned, 0.01);

                // Confirm visibility
                t_ilm_bool visibility_rtn;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetVisibility(layerIDs[j], &visibility_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.visibility,
                            visibility_rtn);

                // Change something
                layer = layers_allocated[index].layerId;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          layers_allocated[index].layerProperties.visibility));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiDoNotSendNotificationsAfterRemoveLayer)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 22 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 42 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 202 * (i + 1);
        layer->layerProperties.origSourceHeight = 244 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        layers_allocated[i].layerProperties.orientation = orientations[i % 4];
        callbackLayerId = layers_allocated[i].layerId;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertCallbackcalled();
        ASSERT_EQ(callbackLayerId,layers_allocated[i].layerId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            callbackLayerId = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Try to modify removed layers orientation
        // Verify no callback is made when change attempted
        ASSERT_EQ(ILM_FAILED,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();      

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                          dimreturned[1]);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_layerGetOrientation(layerIDs[j],
                          &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // Expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multipleRegistrationsLayer)
{
    uint no_surfaces = 8;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 31 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 42 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 198 * (i + 1);
        layer->layerProperties.origSourceHeight = 232 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = { 22 + (i * 5), 42 + (i * 5) };
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Add notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        callbackLayerId = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.orientation = orientations[i % 4];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertCallbackcalled();
        EXPECT_EQ(layers_allocated[i].layerId, callbackLayerId);
    }

    // Remove notification from one layer
    // Change orientation on layer
    // Confirm no callback is made
    // Add callback back, confirm it's back
    for (int i = 0; i < no_layers; i++)
    {
        callbackLayerId = layers_allocated[i].layerId;
        // At the minute this seems to effect the other layer as well - error
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();

        for (int j = i + 1; j < no_layers; j++)
        {
            // Confirm another layer still has call back
            callbackLayerId = layers_allocated[j].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetOrientation(layers_allocated[j].layerId,
                      layers_allocated[j].layerProperties.orientation));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertCallbackcalled(layers_allocated[j].layerId);
            EXPECT_EQ(layers_allocated[j].layerId, callbackLayerId);
        }
    }

    // Return notifications to layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        // add notification
std::cout << "Adding Notification for layer Id: " << layers_allocated[i].layerId << std::endl;
        ilmErrorTypes status = ilm_layerAddNotification(layers_allocated[i].layerId,&LayerCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // Add surfaces to layers check notifications are back
    for (int i = 0; i < no_layers; i++)
    {
        callbackLayerId = layers_allocated[i].layerId;

        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetOrientation(layers_allocated[i].layerId,
                      layers_allocated[i].layerProperties.orientation));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect callback to have been called
            assertCallbackcalled();
            EXPECT_EQ(layers_allocated[i].layerId, callbackLayerId);
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);


                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_layerGetOrientation(layers_allocated[index].layerId,
                          &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Change something
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          layers_allocated[index].layerProperties.visibility));
                ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

                // expect callback to have been called
                assertCallbackcalled();
                EXPECT_EQ(layers_allocated[index].layerId,callbackLayerId);
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiDefaultIsNotToReceiveNotificationsLayer)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 28 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 49 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                  surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 198 * (i + 1);
        layer->layerProperties.origSourceHeight = 232 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = { 22 + (i * 5), 42 + (i * 5) };
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set source rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.sourceX = 3 + (i * 5);
        layers_allocated[i].layerProperties.sourceY = 2 + (i * 5);
        layers_allocated[i].layerProperties.sourceWidth = 50 + (i * 10);
        layers_allocated[i].layerProperties.sourceHeight = 55 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.destX = 5 + (i * 5);
        layers_allocated[i].layerProperties.destY = 3 + (i * 5);
        layers_allocated[i].layerProperties.destWidth = 51 + (i * 10);
        layers_allocated[i].layerProperties.destHeight = 63 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.destX,
                                              layers_allocated[i].layerProperties.destY,
                                              layers_allocated[i].layerProperties.destWidth,
                                              layers_allocated[i].layerProperties.destHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        layer = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.orientation = orientations[i % 4];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set Opacity of layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.opacity = 0.23 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOpacity(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set Visibility
    for (int i = 0; i < no_layers; i++)
    {
        layers_allocated[i].layerProperties.visibility = visibility[no_layers % 2];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            layer = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect no callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // Change position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = { 19 + (i * 5), 37 + (i * 5) };
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                  layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect no callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set source rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.sourceX = 7 + (i * 5);
        layers_allocated[i].layerProperties.sourceY = 5 + (i * 5);
        layers_allocated[i].layerProperties.sourceWidth = 58 + (i * 10);
        layers_allocated[i].layerProperties.sourceHeight = 59 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set destination rectangle and check callback
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer = layers_allocated[i].layerId;

        layers_allocated[i].layerProperties.destX = 8 + (i * 5);
        layers_allocated[i].layerProperties.destY = 2 + (i * 5);
        layers_allocated[i].layerProperties.destWidth = 64 + (i * 10);
        layers_allocated[i].layerProperties.destHeight = 69 + (i * 10);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetDestinationRectangle(layers_allocated[i].layerId,
                                                   layers_allocated[i].layerProperties.destX,
                                                   layers_allocated[i].layerProperties.destY,
                                                   layers_allocated[i].layerProperties.destWidth,
                                                   layers_allocated[i].layerProperties.destHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        // expect callback to have been called
        assertNoCallbackIsCalled();
    }

    // Set Orientation of layers check callbacks
    for (int i = 0; i < no_layers; i++)
    {
        layer = layers_allocated[i].layerId;
        layers_allocated[i].layerProperties.orientation
              = orientations[3 - (i % 4)];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set Opacity of layers
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.opacity = 0.37 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOpacity(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    // Set Visibility
    for (int i = 0; i < no_layers; i++)
    {
        layers_allocated[i].layerProperties.visibility
                  = visibility[1 - (no_layers % 2)];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }


    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId) index = k;
            }

            std::cout << "index set to: " << index << std::endl;

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                callbackLayerId = layers_allocated[index].layerId;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                EXPECT_EQ(layers_allocated[index].layerProperties.destWidth, dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.destHeight, dimreturned[1]);

                // Confirm source rectangle
                ilmLayerProperties layerProperties;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                EXPECT_EQ(layers_allocated[index].layerProperties.sourceX,
                          layerProperties.sourceX);
                EXPECT_EQ(layers_allocated[index].layerProperties.sourceY,
                          layerProperties.sourceY);
                EXPECT_EQ(layers_allocated[index].layerProperties.sourceWidth,
                          layerProperties.sourceWidth);
                EXPECT_EQ(layers_allocated[index].layerProperties.sourceHeight,
                          layerProperties.sourceHeight);

                // Confirm destination rectangle
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfLayer(layers_allocated[index].layerId,
                                                   &layerProperties));
                EXPECT_EQ(layers_allocated[index].layerProperties.destX,
                          layerProperties.destX);
                EXPECT_EQ(layers_allocated[index].layerProperties.destY,
                          layerProperties.destY);
                EXPECT_EQ(layers_allocated[index].layerProperties.destWidth,
                          layerProperties.destWidth);
                EXPECT_EQ(layers_allocated[index].layerProperties.destHeight,
                          layerProperties.destHeight);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_layerGetOrientation(layers_allocated[index].layerId,
                          &orientation_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                          orientation_rtn);

                // Confirm opacity
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layerIDs[j], &returned));
                EXPECT_NEAR(layers_allocated[index].layerProperties.opacity,
                            returned, 0.01);

                // Confirm visibility
                t_ilm_bool visibility_rtn;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetVisibility(layerIDs[j], &visibility_rtn));
                EXPECT_EQ(layers_allocated[index].layerProperties.visibility,
                            visibility_rtn);

                // Change something
                layer = layers_allocated[index].layerId;
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          layers_allocated[index].layerProperties.visibility));

                // expect callback to have been called
                assertNoCallbackIsCalled();
            }
        }
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_layer> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();
}


TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnSurfaceSetPosition)
{
    uint no_surfaces = 4;
    uint no_layers = 4;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 15 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 25 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Add notifications to surface
    for (uint i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceAddNotification(surfaces_allocated[i].returnedSurfaceId,
                  &SurfaceCallbackFunction));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces check call backs
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_pos[2] = {22 + (i * 5), 42 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                  surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        callbackSurfaceId = surfaces_allocated[i].returnedSurfaceId;                
        assertCallbackcalled();
        EXPECT_EQ(callbackSurfaceId,
                  surfaces_allocated[i].returnedSurfaceId);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm parameters are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (IDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceX,
                          posreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceY,
                          posreturned[1]);

            }
        }

        free(IDs);
    }

    surfaces_allocated.clear();
}

TEST_F(IlmMultiNotificationTest, ilm_multiNotifyOnSurfaceSetDimension)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 25 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 15 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        surfaces_allocated.push_back(*surface);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 201 * (i + 1);
        layer->layerProperties.origSourceHeight = 202 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        layers_allocated.push_back(*layer);
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                  surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // assert that we have not been notified
    assertNoCallbackIsCalled();

    // Add notifications to surface
    for (uint i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceAddNotification(surfaces_allocated[i].returnedSurfaceId,
                  &SurfaceCallbackFunction));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimesion again and check call backs
    for (uint i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                  surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        callbackSurfaceId = surfaces_allocated[i].returnedSurfaceId;                
        assertCallbackcalled();
        EXPECT_EQ(callbackSurfaceId,
                  surfaces_allocated[i].returnedSurfaceId);
    }

    // Add surfaces to layers check notifications
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = i * (no_surfaces / no_layers);
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            layer = layers_allocated[i].layerId;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
            // expect no callback to have been called
            assertNoCallbackIsCalled();
        }
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());
        assertNoCallbackIsCalled();
    }

    layers_allocated.clear();

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm parameters are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_surfaces;

            for (int k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId) index = k;
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);
                // Set the dimesnions again and check the callback
                t_ilm_uint surf_dim[2] = {surfaces_allocated[index].surfaceProperties.origSourceWidth,
                                          surfaces_allocated[index].surfaceProperties.origSourceHeight};
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceSetDimension(surfaces_allocated[index].returnedSurfaceId,
                          surf_dim));
                EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());
                callbackSurfaceId = surfaces_allocated[index].returnedSurfaceId;                
                assertCallbackcalled();
                EXPECT_EQ(callbackSurfaceId,
                          surfaces_allocated[index].returnedSurfaceId);

            }
        }
    }

    surfaces_allocated.clear();
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

    bool exit_fail = false;

    std::cout << "Number of parameters: " << argc << std::endl;  

    while ((argc > 1) && (argv[1][0] == '-') && (argv[1][1] == '-'))
    {
        switch (argv[1][2])
        {
            case 'l':
                {
                    std::string layerID(&argv[1][3]);
                    std::istringstream iss(layerID);
                    unsigned int ID;
                    iss >> ID;
                    if (!iss.fail())
                    {
                        TestBase::setStartLayerId(ID);
                        std::cout << "Starting Layer ID set to: " << ID << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid layer start parameter";
                    }
                }
            break;
            case 's':
                {
                    std::string surfaceID(&argv[1][3]);
                    std::istringstream iss(surfaceID);
                    unsigned short int ID;
                    iss >> ID;
                    if (!iss.fail())
                    {
                        TestBase::setStartSurfaceId(ID);
                        std::cout << "Starting Surface ID set to: " << ID << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid surface start parameter";
                    }
                }
            break;
            case 'n':
                {
                    std::string layers(&argv[1][3]);
                    std::istringstream iss(layers);
                    unsigned short int maxNumberLayers;
                    iss >> maxNumberLayers;
                    if (!iss.fail())
                    {
                        TestBase::setMaxLayerIds(maxNumberLayers);
                        std::cout << "Maximum number of contiguous layers set to: " << maxNumberLayers << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid maximum layer parameter";
                    }
                }
            break;
            case 't':
                {
                    std::string layers(&argv[1][3]);
                    std::istringstream iss(layers);
                    unsigned short int maxNumberSurfaces;
                    iss >> maxNumberSurfaces;
                    if (!iss.fail())
                    {
                        TestBase::setMaxSurfaceIds(maxNumberSurfaces);
                        std::cout << "Maximum number of contiguous surfaces set to: " << maxNumberSurfaces << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid maximum surface parameter";
                    }
                }
            break;

            default:
            {
                 EXPECT_TRUE(false) << "Unknown parameter specified: " << argv[1][3] << std::endl;
                 EXPECT_TRUE(false) << "Options: --l<starting layer id>" << std::endl
                                    << "         --s<starting surface id>" << std::endl
                                    << "         --n<number of layers to be used>" << std::endl
                                    << "         --t<number of surfaces to be used>" << std::endl;
                 exit_fail = true;
            }
        }

        ++argv;
        --argc;
    }

    if (exit_fail)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return RUN_ALL_TESTS();
    }
}
