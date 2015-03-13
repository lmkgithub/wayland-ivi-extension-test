/***************************************************************************
 *
 * Copyright 2015 Codethink Ltd
 * Copyright 2010-2014 BMW Car IT GmbH
 * Copyright (C) 2012 DENSO CORPORATION and Robert Bosch Car Multimedia Gmbh
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

class IlmNullPointerTest : public TestBase, public ::testing::Test {

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

    IlmNullPointerTest(){}
    ~IlmNullPointerTest(){}

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
t_ilm_layer IlmNullPointerTest::callbackLayerId;
t_ilm_surface IlmNullPointerTest::callbackSurfaceId;
struct ilmLayerProperties IlmNullPointerTest::LayerProperties;
unsigned int IlmNullPointerTest::mask;
t_ilm_surface IlmNullPointerTest::surface;
ilmSurfaceProperties IlmNullPointerTest::SurfaceProperties;

TEST_F(IlmNullPointerTest, ilm_testSurfaceCreateDimensionNullPointer)
{
    uint no_surfaces = 4;

    // Try to create surfaces using null ptr for id ref.
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 150 + ( i * 10 );
        surface->surfaceProperties.origSourceHeight = 250 + ( i + 10 );
        surfaces_allocated.push_back(*surface);

        ASSERT_EQ(ILM_FAILED, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     NULL));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Now try to create surfaces with valid id ref.
    for (int i = 0; i < surfaces_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                     surfaces_allocated[i].surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_allocated[i].returnedSurfaceId)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces using null pointer
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_FAILED, ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId, NULL));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces using valid pointer
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check parameters
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dim_rtn[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                  dim_rtn));

        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth,
                  dim_rtn[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight,
                  dim_rtn[1]);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
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
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint surf_pos[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_surfaceSetVisibility(surfaces_allocated[index].returnedSurfaceId,
                          ILM_TRUE));

                // expect callback to have been called
                assertNoCallbackIsCalled();
            }
        }

        surfaceIDs.clear();
    }

    surfaces_allocated.clear();
}

TEST_F(IlmNullPointerTest, ilm_testLayerCreateDimensionNullPointer)
{
    uint no_layers = 4;

    // Try to create layers using null ptr for id ref.
    for (int i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 150 + ( i * 10 );
        layer->layerProperties.origSourceHeight = 250 + ( i + 10 );
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_FAILED,
                  ilm_layerCreateWithDimension(NULL,
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Now try to create layers with valid id ref.
    for (int i = 0; i < layers_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers_allocated[i].layerId),
                                               layers_allocated[i].layerProperties.origSourceWidth,
                                               layers_allocated[i].layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of layers using null pointer
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_FAILED,
                  ilm_layerSetDimension(layers_allocated[i].layerId, NULL));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of layers using valid pointer
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {layers_allocated[i].layerProperties.origSourceWidth,
                                  layers_allocated[i].layerProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetDimension(layers_allocated[i].layerId,
                                        surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check parameters
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint dim_rtn[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetDimension(layers_allocated[i].layerId,
                  dim_rtn));

        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceWidth, dim_rtn[0]);
        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceHeight, dim_rtn[1]);
    }

    // Loop through layers and remove
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining layers and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            uint index = no_layers;

            for (int k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != no_layers)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint surf_pos[2] = {0, 0};
                ilmLayerProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                          dimreturned[1]);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_layerSetVisibility(layers_allocated[index].layerId,
                          ILM_TRUE));

                // expect callback to have been called
                assertNoCallbackIsCalled();
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

TEST_F(IlmNullPointerTest, ilm_testSurfaceSetPositionNullPointer)
{
    uint no_surfaces = 4;

    // Try to create surfaces using null ptr for id ref.
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 150 + ( i * 10 );
        surface->surfaceProperties.origSourceHeight = 250 + ( i + 10 );
        surfaces_allocated.push_back(*surface);

        ASSERT_EQ(ILM_FAILED, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     NULL));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Now try to create surfaces with valid id ref.
    for (int i = 0; i < surfaces_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                     surfaces_allocated[i].surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_allocated[i].returnedSurfaceId)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces using null pointer
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_FAILED,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          NULL));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces using valid pointer
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check parameters
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dim_rtn[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          dim_rtn));

        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth,
                  dim_rtn[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight,
                  dim_rtn[1]);
    }

    // Try to set position of surfaces using null pointer
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_FAILED,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         NULL));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of surfaces using valid pointer
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {10, 15};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check position of surfaces set correctly
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(surf_pos[0], surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(surf_pos[1], surfaces_allocated[i].surfaceProperties.sourceY);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;
        std::vector<t_ilm_surface> surfaceIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceRemoveNotification(surfaces_allocated[i].returnedSurfaceId));
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceRemove(surfaces_allocated[i].returnedSurfaceId));
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
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint surf_pos[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                // Change something that has been pre-set and check callback
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_surfaceSetVisibility(surfaces_allocated[index].returnedSurfaceId,
                          ILM_TRUE));

                // expect callback to have been called
                assertNoCallbackIsCalled();
            }
        }

        surfaceIDs.clear();
    }

    surfaces_allocated.clear();
}
