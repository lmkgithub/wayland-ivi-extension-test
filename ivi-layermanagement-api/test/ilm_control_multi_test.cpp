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

class IlmCommandMultiTest : public TestBase, public ::testing::Test {

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

    IlmCommandMultiTest(){}
    ~IlmCommandMultiTest(){}

    void SetUp()
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_initWithNativedisplay((t_ilm_nativedisplay)wlDisplay));
        mask = static_cast<t_ilm_notification_mask>(0);
        surface = INVALID_ID;
        timesCalled=0;

        callbackLayerId = INVALID_ID;
        callbackSurfaceId = INVALID_ID;
    }

    void TearDown()
    {
        // set default values
        callbackLayerId = INVALID_ID;
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

        EXPECT_TRUE(callbackLayerId == INVALID_ID || callbackLayerId == layer);
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

//        EXPECT_TRUE(callbackSurfaceId == INVALID_ID || callbackSurfaceId == surface);
        callbackSurfaceId = surface;
        mask |= (unsigned)m;
        timesCalled++;

        pthread_cond_signal( &waiterVariable );
    }
};

// Pointers where to put received values for current Test
t_ilm_layer IlmCommandMultiTest::callbackLayerId;
t_ilm_surface IlmCommandMultiTest::callbackSurfaceId;
struct ilmLayerProperties IlmCommandMultiTest::LayerProperties;
unsigned int IlmCommandMultiTest::mask;
t_ilm_surface IlmCommandMultiTest::surface;
ilmSurfaceProperties IlmCommandMultiTest::SurfaceProperties;

/*****************************************************************************
 * This test checks the setting and getting of surface dimensions when       *
 * multiple surfaces are used and they are associated with layers.           *
 * Multiple surfaces are created and their dimensions are set. Multiple      *
 * layers are created and then the surfaces are added to the layers.         *
 * The dimensions for the surfaces are checked to make sure they are         *
 * consistent with those originally set.                                     *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions as originally set.          *
 *****************************************************************************/

TEST_F(IlmCommandMultiTest, multi_SetGetSurfaceDimension) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
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
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Get dimensions for surfaces and compare against those originally set
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dimreturned[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                                              dimreturned));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth,
                  dimreturned[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight,
                  dimreturned[1]);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                callbackSurfaceId = surfaces_allocated[index].returnedSurfaceId;

                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of layer dimensions when         *
 * multiple layers are used and surfaces are associated with them.           *
 * Multiple surfaces are created and their dimensions are set. Multiple      *
 * layers are created and then the surfaces are added to the layers.         *
 * The dimensions for the layers are checked to make sure they are           *
 * consistent with those originally set.                                     *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions as originally set.          *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetLayerDimension) {
    uint no_surfaces = 16;
    uint no_layers = 4;

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
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
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                 j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Get dimensions for layers and compare against those originally set
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint dimreturned[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetDimension(layers_allocated[i].layerId,
                                        dimreturned));
        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceWidth,
                  dimreturned[0]);
        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceHeight,
                  dimreturned[1]);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                ilmSurfaceProperties returnValue;
                callbackSurfaceId = surfaces_allocated[index].returnedSurfaceId;

                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth, dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight, dimreturned[1]);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of surface positions when        *
 * multiple layers are used and surfaces are associated with them. Multiple  *
 * surfaces are created and their dimensions and positions are set. Multiple *
 * layers are created (dimensions and positions set) and then the surfaces   *
 * are added to the layers.                                                  *
 * The positions for the surfaces and layers are checked to make sure they   *
 * are consistent with those originally set.                                 *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions as originally set.          *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetSurfacePosition)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
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
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {15 + (i * 5), 25 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {25 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Check position of layers set correctly
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
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

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaceIDs[j],
                          posreturned));
                EXPECT_EQ(posreturned[0],
                          surfaces_allocated[index].surfaceProperties.sourceX);
                EXPECT_EQ(posreturned[1],
                          surfaces_allocated[index].surfaceProperties.sourceY);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint layer_pos[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                   layer_pos));
                    ASSERT_EQ(layer_pos[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(layer_pos[1],
                              layers_allocated[index].layerProperties.sourceY);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of layer positions when          *
 * multiples are used and surfaces are associated with them. Multiple        *
 * surfaces are created and their dimensions and positions are set. Multiple *
 * layers are created (dimensions and positions set) and then the surfaces   *
 * are added to the layers.                                                  *
 * The positions for the surfaces and layers are checked to make sure they   *
 * are consistent with those originally set.                                 *
 * The position of layers are changed and then re-checked.                   *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions as originally set.          *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetLayerPosition) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
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
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {15 + (i * 5), 25 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {25 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }

    // Change position of layers again
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {32 + (i * 15), 25 + (i * 15)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Check position of layers still set correctly
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }

    // Check position of surfaces still set correctly
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

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaceIDs[j],
                          posreturned));
                EXPECT_EQ(posreturned[0],
                          surfaces_allocated[index].surfaceProperties.sourceX);
                EXPECT_EQ(posreturned[1],
                          surfaces_allocated[index].surfaceProperties.sourceY);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint layer_pos[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                   layer_pos));
                    ASSERT_EQ(layer_pos[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(layer_pos[1],
                              layers_allocated[index].layerProperties.sourceY);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of surface orientations when     *
 * multiples are used and they are associated with layers. Multiple          *
 * surfaces are created and their dimensions and orientations are set.       *
 * Layers are created (dimensions set).                                      *
 * The positions for the surfaces and layers are set and checked.            *
 * The surfaces are added to the layers. The surface orientations are re-    *
 * checked to confirm that they are still valid.                             *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions and orientations (surfaces  *
 * only) as originally set.                                                  *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetSurfaceOrientation) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
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
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Set Orientation of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.orientation
            = orientations[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  surfaces_allocated[i].surfaceProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Orientation of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                            &orientation_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.orientation,
                  orientation_rtn);
    }

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {15 + (i * 5), 25 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {25 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Confirm Orientation of surfaces still valid
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                            &orientation_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.orientation,
                  orientation_rtn);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                e_ilmOrientation orientation_rtn;

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaceIDs[j],
                          posreturned));
                EXPECT_EQ(posreturned[0],
                          surfaces_allocated[index].surfaceProperties.sourceX);
                EXPECT_EQ(posreturned[1],
                          surfaces_allocated[index].surfaceProperties.sourceY);

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetOrientation(surfaceIDs[j],
                                                    &orientation_rtn));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.orientation,
                          orientation_rtn);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint layer_pos[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                   layer_pos));
                    ASSERT_EQ(layer_pos[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(layer_pos[1],
                              layers_allocated[index].layerProperties.sourceY);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of layer orientations when       *
 * multiples are used and they are surfaces are associated with them.        *
 * Multiple surfaces are created and their dimensions  are set.              *
 * Layers are created with dimensions set.                                   *
 * The positions for the surfaces and layers are set and checked.            *
 * The orientations of layers are set and confirmed.                         *
 * The surfaces are added to the layers. The layer orientations are re-      *
 * checked to confirm that they are still valid.                             *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions, positions and orientations *
 * (layers only) are as originally set.                                      *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetLayerOrientation) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientations[4] = {ILM_ZERO, ILM_NINETY,
                                        ILM_ONEHUNDREDEIGHTY,
                                        ILM_TWOHUNDREDSEVENTY};

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 17 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 23 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {15 + (i * 5), 25 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {25 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }


    // Set Orientation of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.orientation
            = orientations[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Orientation of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                                          &orientation_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.orientation,
                  orientation_rtn);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Confirm Orientation of layers still valid
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                                          &orientation_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.orientation,
                  orientation_rtn);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                e_ilmOrientation orientation_rtn;

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaceIDs[j],
                          posreturned));
                EXPECT_EQ(posreturned[0],
                          surfaces_allocated[index].surfaceProperties.sourceX);
                EXPECT_EQ(posreturned[1],
                          surfaces_allocated[index].surfaceProperties.sourceY);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint layer_pos[2] = {0, 0};
                    e_ilmOrientation orientation_rtn;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                   layer_pos));
                    ASSERT_EQ(layer_pos[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(layer_pos[1],
                              layers_allocated[index].layerProperties.sourceY);


                    // Confirm orientation
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetOrientation(layerIDs[j],
                                                      &orientation_rtn));
                    EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                              orientation_rtn);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of multiple layer and surface    *
 * opacities (this essentially combines two tests for brevity).              *
 * Multiple surfaces are created and their dimension and opacities are set.  *
 * Layers are created with dimensions set.                                   *
 * The positions for the surfaces and layers are set and checked.            *
 * The opacities of layers are set and confirmed.                            *
 * The surfaces are added to the layers. The layer & surface opacities are   *
 * re-checked to confirm that they are still valid.                          *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions, positions and opacities    *
 * are as originally set.                                                    *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetLayerSurfaceOpacity) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 17 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 23 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Set Opacity of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.opacity = 0.05 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOpacity(surfaces_allocated[i].returnedSurfaceId,
                                        surfaces_allocated[i].surfaceProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Opacity of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_float returned;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOpacity(surfaces_allocated[i].returnedSurfaceId,
                  &returned));
        EXPECT_NEAR(surfaces_allocated[i].surfaceProperties.opacity,
                    returned, 0.01);
    }

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {15 + (i * 5), 25 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {25 + (i * 5), 40 + (i * 5)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }


    // Set Opacity of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.opacity = 0.05 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOpacity(layers_allocated[i].layerId,
                                      layers_allocated[i].layerProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Opacity of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_float returned;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOpacity(layers_allocated[i].layerId,
                  &returned));
        EXPECT_NEAR(layers_allocated[i].layerProperties.opacity,
                    returned, 0.01);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Confirm Opacity of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_float returned;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOpacity(layers_allocated[i].layerId,
                  &returned));
        EXPECT_NEAR(layers_allocated[i].layerProperties.opacity,
                    returned, 0.01);
    }

    // Confirm Opacity of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_float returned;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOpacity(surfaces_allocated[i].returnedSurfaceId,
                  &returned));
        EXPECT_NEAR(surfaces_allocated[i].surfaceProperties.opacity,
                    returned, 0.01);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                e_ilmOrientation orientation_rtn;

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                // Confirm opacity
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetOpacity(surfaceIDs[j],
                          &returned));
                EXPECT_NEAR(surfaces_allocated[index].surfaceProperties.opacity,
                            returned, 0.01);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint layer_pos[2] = {0, 0};
                    e_ilmOrientation orientation_rtn;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                   layer_pos));
                    ASSERT_EQ(layer_pos[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(layer_pos[1],
                              layers_allocated[index].layerProperties.sourceY);


                    // Confirm opacity
                    t_ilm_float returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetOpacity(layerIDs[j],
                                                  &returned));
                    EXPECT_NEAR(layers_allocated[index].layerProperties.opacity,
                                returned, 0.01);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of surface visibilities when     *
 * multiples are used and they are associated with layers. Multiple          *
 * surfaces are created and their dimensions and visibilities are set.       *
 * Layers are created (dimensions and visibilities set).                     *
 * The positions for the surfaces are set.                                   *
 * The surfaces are added to the layers. The surface and layer visibilities  *
 * are rechecked to confirm that they are still valid.                       *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions and visibilities (surfaces  *
 * only) as originally set.                                                  *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetGetSurfaceVisibility) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 17 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 23 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Set Visibility
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.visibility
            = visibility[no_layers % 2];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetVisibility(surfaces_allocated[i].returnedSurfaceId,
                  surfaces_allocated[i].surfaceProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Visibility
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetVisibility(surfaces_allocated[i].returnedSurfaceId,
                  &visibility_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.visibility,
                  visibility_rtn);
    }

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.visibility
            = visibility[no_layers % 2];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Visibility
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetVisibility(layers_allocated[i].layerId,
                  &visibility_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                  visibility_rtn);
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {15 + (i * 5), 25 + (i * 5)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
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

    // Confirm Visibility is still valid
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetVisibility(surfaces_allocated[i].returnedSurfaceId,
                  &visibility_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.visibility,
                  visibility_rtn);
    }

    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetVisibility(layers_allocated[i].layerId,
                  &visibility_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                  visibility_rtn);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                t_ilm_bool visibility_rtn;

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                // Confirm visibility
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetVisibility(surfaces_allocated[index].returnedSurfaceId,
                                                   &visibility_rtn));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.visibility,
                          visibility_rtn);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_bool visibility_rtn;

                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm visibility
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetVisibility(layers_allocated[index].layerId,
                                                     &visibility_rtn));
                    EXPECT_EQ(layers_allocated[index].layerProperties.visibility,
                              visibility_rtn);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of surface source rectangle when *
 * multiples are used and they are associated with layers. Multiple          *
 * surfaces are created and their dimensions and source rectangles are set.  *
 * Layers are created (dimensions and visibilities set).                     *
 * The surfaces are added to the layers. The surface source rectangles are   *
 * rechecked to confirm that they are still valid.                           *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions and visibilities (surfaces  *
 * only) as originally set.                                                  *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetSurfaceSourceRectangle) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 17 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 23 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Set source rectangle of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.sourceX = 82 + (i * 10);
        surfaces_allocated[i].surfaceProperties.sourceY = 6238 + (i * 3);
        surfaces_allocated[i].surfaceProperties.sourceWidth = 618 + (i * 7);
        surfaces_allocated[i].surfaceProperties.sourceHeight = 3 + (i * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetSourceRectangle(surfaces_allocated[i].returnedSurfaceId,
                                                surfaces_allocated[i].surfaceProperties.sourceX,
                                                surfaces_allocated[i].surfaceProperties.sourceY,
                                                surfaces_allocated[i].surfaceProperties.sourceWidth,
                                                surfaces_allocated[i].surfaceProperties.sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm surface source rectangle
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceX,
                  surfaceProperties.sourceX);
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceY,
                  surfaceProperties.sourceY);
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceWidth,
                  surfaceProperties.sourceWidth);
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceHeight,
                  surfaceProperties.sourceHeight);
    }

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.visibility
            = visibility[no_layers % 2];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Visibility
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetVisibility(layers_allocated[i].layerId,
                  &visibility_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                  visibility_rtn);
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
                  j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }
    }


    // Confirm surface source rectangle still correct
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceX,
                  surfaceProperties.sourceX);
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceY,
                  surfaceProperties.sourceY);
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceWidth,
                  surfaceProperties.sourceWidth);
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceHeight,
                  surfaceProperties.sourceHeight);
    }

    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetVisibility(layers_allocated[i].layerId,
                  &visibility_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                  visibility_rtn);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                t_ilm_bool visibility_rtn;
                ilmSurfaceProperties surfaceProperties;

                // Check dimensions
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);


                // Check rectangle
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfSurface(surfaces_allocated[index].returnedSurfaceId,
                                                     &surfaceProperties));
                ASSERT_EQ(surfaces_allocated[index].surfaceProperties.sourceX,
                          surfaceProperties.sourceX);
                ASSERT_EQ(surfaces_allocated[index].surfaceProperties.sourceY,
                          surfaceProperties.sourceY);
                ASSERT_EQ(surfaces_allocated[index].surfaceProperties.sourceWidth,
                          surfaceProperties.sourceWidth);
                ASSERT_EQ(surfaces_allocated[index].surfaceProperties.sourceHeight,
                          surfaceProperties.sourceHeight);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_bool visibility_rtn;

                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the setting and getting of layer source rectangle when   *
 * multiples are used and they are associated with surfaces. Multiple        *
 * surfaces are created and their dimensions are set.                        *
 * Layers are created (dimensions and visibilities set).                     *
 * The source rectangles of layers are set.                                  *
 * The surfaces are added to the layers. The layer source rectangles are     *
 * rechecked to confirm that they are still valid.                           *
 * The visibilities of layers are rechecked.                                 *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions and source rectangles       *
 * (layers only) as originally set.                                          *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, multi_SetLayerSourceRectangle) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 17 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 23 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.visibility
            = visibility[no_layers % 2];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[i].layerId,
                  layers_allocated[i].layerProperties.visibility));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm Visibility
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetVisibility(layers_allocated[i].layerId,
                  &visibility_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                  visibility_rtn);
    }

    // Set source rectangle of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.sourceX = 82 + (i * 10);
        layers_allocated[i].layerProperties.sourceY = 6238 + (i * 3);
        layers_allocated[i].layerProperties.sourceWidth = 618 + (i * 7);
        layers_allocated[i].layerProperties.sourceHeight = 3 + (i * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.sourceX,
                                              layers_allocated[i].layerProperties.sourceY,
                                              layers_allocated[i].layerProperties.sourceWidth,
                                              layers_allocated[i].layerProperties.sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm layers source rectangle
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmLayerProperties layerProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers_allocated[i].layerId,
                                           &layerProperties));
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceX,
                  layerProperties.sourceX);
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceY,
                  layerProperties.sourceY);
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceWidth,
                  layerProperties.sourceWidth);
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceHeight,
                  layerProperties.sourceHeight);
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers check notifications
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = i * (surfaces_allocated.size() / layers_allocated.size());
                  j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
                  j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                      surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }
    }


    // Confirm layer source rectangle still correct
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmLayerProperties layerProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers_allocated[i].layerId,
                                           &layerProperties));
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceX,
                  layerProperties.sourceX);
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceY,
                  layerProperties.sourceY);
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceWidth,
                  layerProperties.sourceWidth);
        ASSERT_EQ(layers_allocated[i].layerProperties.sourceHeight,
                  layerProperties.sourceHeight);
    }

    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm visibility
        t_ilm_bool visibility_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetVisibility(layers_allocated[i].layerId,
                  &visibility_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                  visibility_rtn);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                t_ilm_bool visibility_rtn;

                // Check dimensions
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_bool visibility_rtn;
                    ilmLayerProperties layerProperties;

                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Check rectangle
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
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the screenshotting of both multiple layers and surfaces. *
 * Multiple surfaces are created and their dimensions and positions are set. *
 * Layers are created (dimensions and positions are set).                    *
 * Screen shots are taken of every single layer.                             *
 * Screen shots are taken of every single surface.                           *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions and positions are as        *
 * originally set.                                                           *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, ilm_multiTakeSurfaceLayerScreenshots) {

    uint no_surfaces = 4;
    uint no_layers = 4;

    const char* outputFile = "/tmp/test.bmp";

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 17 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 23 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {18 + (i * 7), 21 + (i * 4)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 200 * (i + 1);
        layer->layerProperties.origSourceHeight = 240 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers using valid pointer
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {12 + (i * 2), 23 + (i * 6)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop round and take layer shots
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // make sure the file is not there before
        FILE* f = fopen(outputFile, "r");
        if (f!=NULL){
            fclose(f);
            int result = remove(outputFile);
            ASSERT_EQ(0, result);
        }

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_takeLayerScreenshot(outputFile,
                                          layers_allocated[i].layerId));

        sleep(1);
        f = fopen(outputFile, "r");
        ASSERT_TRUE(f!=NULL);
        fclose(f);
        remove(outputFile);
    }

    // Loop round and take surface shots
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // make sure the file is not there before
        FILE* f = fopen(outputFile, "r");
        if (f!=NULL){
            fclose(f);
            int result = remove(outputFile);
            ASSERT_EQ(0, result);
        }

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_takeSurfaceScreenshot(outputFile,
                                            surfaces_allocated[i].returnedSurfaceId));

        sleep(1);
        f = fopen(outputFile, "r");
        ASSERT_TRUE(f!=NULL);
        fclose(f);
        remove(outputFile);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};

                // Check dimensions
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                // Confirm position
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetPosition(surfaces_allocated[index].returnedSurfaceId,
                                                 posreturned));
                ASSERT_EQ(posreturned[0], surfaces_allocated[index].surfaceProperties.sourceX);
                ASSERT_EQ(posreturned[1], surfaces_allocated[index].surfaceProperties.sourceY);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint posreturned[2] = {0, 0};

                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                     posreturned));
                    ASSERT_EQ(posreturned[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(posreturned[1],
                              layers_allocated[index].layerProperties.sourceY);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

/*****************************************************************************
 * This test checks the creation and checking of different types of pixel-   *
 * format surfaces.                                                          *
 * Multiple surfaces are created with different pixel formats.               *
 * The dimensions of the surfaces are set.                                   *
 * Layers are created (dimensions are set).                                  *
 * The positions of layers are set.                                          *
 * The positions of surfaces are set.                                        *
 * The pixel formats of every surface are checked.                           *
 * Layers and surfaces are removed one by one, checking on each occasion     *
 * that those remaining have the same dimensions and positions are as        *
 * originally set.                                                           *
 *****************************************************************************/
TEST_F(IlmCommandMultiTest, ilm_multiLayerSurfaceGetPixelformat) {

    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmPixelFormat formats[8] = {ILM_PIXELFORMAT_R_8,
                                   ILM_PIXELFORMAT_RGB_888,
                                   ILM_PIXELFORMAT_RGBA_8888,
                                   ILM_PIXELFORMAT_RGB_565,
                                   ILM_PIXELFORMAT_RGBA_5551,
                                   ILM_PIXELFORMAT_RGBA_6661,
                                   ILM_PIXELFORMAT_RGBA_4444,
                                   ILM_PIXEL_FORMAT_UNKNOWN};

    // Create surfaces.
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 12 * (i + 1);
        surface->surfaceProperties.origSourceHeight = 27 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     formats[i % 8],
                                     &surface->returnedSurfaceId));
        surfaces_allocated.push_back(*surface);
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

    // Try to create layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 206 * (i + 1);
        layer->layerProperties.origSourceHeight = 232 * (i + 1);
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {18 + (i * 7), 21 + (i * 4)};
        surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
        surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surf_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Try to set position of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {28 + (i * 4), 44 + (i * 15)};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check pixel format of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmPixelFormat pf;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetPixelformat(surfaces_allocated[i].returnedSurfaceId, &pf));
        ASSERT_EQ(formats[i % 8], pf);
    }

    t_ilm_uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_surfaces;

            for (uint k = 0; k < surfaces_allocated.size(); k++)
            {
                if (surfaceIDs[j] == surfaces_allocated[k].returnedSurfaceId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};

                // Check dimensions
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);

                // Confirm position
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetPosition(surfaces_allocated[index].returnedSurfaceId,
                                                 posreturned));
                ASSERT_EQ(posreturned[0], surfaces_allocated[index].surfaceProperties.sourceX);
                ASSERT_EQ(posreturned[1], surfaces_allocated[index].surfaceProperties.sourceY);
            }
        }
    }

    surfaces_allocated.clear();

    uint total_layers = layers_allocated.size();

    // remove the layers
    for (uint i = 0; i < total_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        std::vector<t_ilm_layer> layerIDs;

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerRemoveNotification(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers_allocated[i].layerId));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
        layerIDs.assign(IDs, IDs + length);
        free(IDs);

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (uint j = 0; j < length; j++)
        {

            uint index = total_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != total_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint k = 0; k < length; k++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    t_ilm_uint posreturned[2] = {0, 0};

                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Confirm position
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layers_allocated[index].layerId,
                                                     posreturned));
                    ASSERT_EQ(posreturned[0],
                              layers_allocated[index].layerProperties.sourceX);
                    ASSERT_EQ(posreturned[1],
                              layers_allocated[index].layerProperties.sourceY);
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}