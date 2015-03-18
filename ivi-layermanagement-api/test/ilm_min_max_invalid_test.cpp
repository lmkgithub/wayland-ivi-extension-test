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
#include <limits>


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

class IlmMinMaxInvalidTest : public TestBase, public ::testing::Test {

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

    IlmMinMaxInvalidTest(){}
    ~IlmMinMaxInvalidTest(){}

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

        EXPECT_TRUE(callbackLayerId == (unsigned)INVALID_ID || callbackLayerId == layer);
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

//        EXPECT_TRUE(callbackSurfaceId == (unsigned)INVALID_ID || callbackSurfaceId == surface);
        callbackSurfaceId = surface;
        mask |= (unsigned)m;
        timesCalled++;

        pthread_cond_signal( &waiterVariable );
    }
};

// Pointers where to put received values for current Test
t_ilm_layer IlmMinMaxInvalidTest::callbackLayerId;
t_ilm_surface IlmMinMaxInvalidTest::callbackSurfaceId;
struct ilmLayerProperties IlmMinMaxInvalidTest::LayerProperties;
unsigned int IlmMinMaxInvalidTest::mask;
t_ilm_surface IlmMinMaxInvalidTest::surface;
ilmSurfaceProperties IlmMinMaxInvalidTest::SurfaceProperties;

TEST_F(IlmMinMaxInvalidTest, ilm_maxminSetGetSurfaceLayerDimension)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def* surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 0;
        surface->surfaceProperties.origSourceHeight = 0;
        surfaces_allocated.push_back(*surface);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
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

    // Get dimensions for surfaces and compare against those originally set
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dim_rtn[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                  dim_rtn));

        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth, dim_rtn[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight, dim_rtn[1]);
    }

    // Create layers
    for (uint i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 0;
        layer->layerProperties.origSourceHeight = 0;
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for layers and compare against those originally set
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint dim_rtn[2] = {1, 1};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetDimension(layers_allocated[i].layerId,
                  dim_rtn));

        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceWidth,
                  dim_rtn[0]);
        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceHeight,
                  dim_rtn[1]);
    }

    // Change dimensions of layers to maximum
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_dim[2];
        layers_allocated[i].layerProperties.origSourceWidth
            = std::numeric_limits<t_ilm_uint>::max();
        layers_allocated[i].layerProperties.origSourceHeight
            = std::numeric_limits<t_ilm_uint>::max();
        layer_dim[0] = layers_allocated[i].layerProperties.origSourceWidth;
        layer_dim[1] = layers_allocated[i].layerProperties.origSourceHeight;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetDimension(layers_allocated[i].layerId, layer_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for layers and compare against modified values
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layers_allocated[i].layerId, dimreturned));
        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceWidth, dimreturned[0]);
        EXPECT_EQ(layers_allocated[i].layerProperties.origSourceHeight, dimreturned[1]);
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2];
        surf_dim[0] = std::numeric_limits<t_ilm_uint>::max();
        surf_dim[1] = std::numeric_limits<t_ilm_uint>::max();
        surfaces_allocated[i].surfaceProperties.origSourceWidth = surf_dim[0];
        surfaces_allocated[i].surfaceProperties.origSourceHeight = surf_dim[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for surfaces and compare
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dim_rtn[2] = {0, 0};

        // Confirm dimension
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                  dim_rtn));

        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth, dim_rtn[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight, dim_rtn[1]);
    }

    // Check and clear any previous surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId,
                                           &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = (i * (surfaces_allocated.size() / layers_allocated.size()));
             j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                                          surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }
    }

    // Get dimensions for surfaces and compare against those originally set
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                                         dimreturned));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth,
                  dimreturned[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight,
                  dimreturned[1]);
    }

    uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetDimension(surfaceIDs[j], dimreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceWidth,
                          dimreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.origSourceHeight,
                          dimreturned[1]);
            }
        }

        surfaceIDs.clear();
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
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    // Change something that has been pre-set and check callback
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerSetVisibility(layers_allocated[index].layerId,
                              ILM_TRUE));

                    // expect callback to not have been called
                    assertNoCallbackIsCalled();
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
};

TEST_F(IlmMinMaxInvalidTest, ilm_maxminSetGetSurfaceLayerPosition)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    // Create surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 150 + ( i * 10 );
        surface->surfaceProperties.origSourceHeight = 250 + ( i + 10 );
        surfaces_allocated.push_back(*surface);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (uint i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = 160 + ( i * 10 );
        layer->layerProperties.origSourceHeight = 170 + ( i * 10 );
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces to maximum
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {std::numeric_limits<t_ilm_uint>::max(),
                                  std::numeric_limits<t_ilm_uint>::max()};
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
        t_ilm_uint surface_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surface_pos));
        ASSERT_EQ(surface_pos[0],
                  surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(surface_pos[1],
                  surfaces_allocated[i].surfaceProperties.sourceY);
    }

    // Set position of layers to maximum
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {std::numeric_limits<t_ilm_uint>::max(),
                                   std::numeric_limits<t_ilm_uint>::max()};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
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

    uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetPosition(surfaceIDs[j], posreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceX,
                          posreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceY,
                          posreturned[1]);
            }
        }

        surfaceIDs.clear();
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
                    t_ilm_uint posreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layerIDs[j], posreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.sourceX,
                              posreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.sourceY,
                              posreturned[1]);

                    // Change something that has been pre-set and check callback
                    ASSERT_EQ(ILM_SUCCESS,
                              ilm_layerSetVisibility(layers_allocated[index].layerId,
                              ILM_TRUE));

                    // expect callback to have been called
                    assertNoCallbackIsCalled();
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}

TEST_F(IlmMinMaxInvalidTest, ilm_maxminSetGetSurfaceLayerOrientation)
{
    uint no_surfaces = 4;
    uint no_layers = 2;

    e_ilmOrientation orientation[4] = {ILM_ZERO, ILM_NINETY,
                                       ILM_ONEHUNDREDEIGHTY,
                                       ILM_TWOHUNDREDSEVENTY};

    // Create surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        surface_def * surface = new surface_def;
        surface->requestedSurfaceId = getSurface();
        surface->returnedSurfaceId = surface->requestedSurfaceId;
        surface->surfaceProperties.origSourceWidth = 150 + ( i * 10 );
        surface->surfaceProperties.origSourceHeight = 250 + ( i + 10 );
        surfaces_allocated.push_back(*surface);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surface->surfaceProperties.origSourceWidth,
                                     surface->surfaceProperties.origSourceHeight,
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surface->returnedSurfaceId)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Orientations of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.orientation = orientation[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                            surfaces_allocated[i].surfaceProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &orientation_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.orientation,
                  orientation_rtn);
    }

    // Create layers
    for (uint i = 0; i < no_layers; i++)
    {
        layer_def * layer = new layer_def;
        layer->layerId = getLayer();
        layer->layerProperties.origSourceWidth = std::numeric_limits<t_ilm_uint>::max();
        layer->layerProperties.origSourceHeight = std::numeric_limits<t_ilm_uint>::max();
        layers_allocated.push_back(*layer);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layer->layerId),
                                               layer->layerProperties.origSourceWidth,
                                               layer->layerProperties.origSourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add surfaces to layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        for (uint j = (i * (surfaces_allocated.size() / layers_allocated.size()));
             j < ((i + 1) * (surfaces_allocated.size() / layers_allocated.size()));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers_allocated[i].layerId,
                                          surfaces_allocated[j].returnedSurfaceId));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }

    }

    // Set Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.orientation = orientation[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                                          layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &orientation_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.orientation,
                  orientation_rtn);
    }

    // Set position of surfaces to maximum
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {std::numeric_limits<t_ilm_uint>::max(),
                                  std::numeric_limits<t_ilm_uint>::max()};
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
        t_ilm_uint surface_pos[2] = {0, 0};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surface_pos));
        ASSERT_EQ(surface_pos[0],
                  surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(surface_pos[1],
                  surfaces_allocated[i].surfaceProperties.sourceY);
    }

    // Set position of layers to maximum
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {std::numeric_limits<t_ilm_uint>::max(),
                                  std::numeric_limits<t_ilm_uint>::max()};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
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

    // Change things again

    // Set Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.orientation = orientation[ 3 - (i % 4) ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                                          layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &orientation_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.orientation,
                  orientation_rtn);
    }

    // Set position of surfaces to minimum
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {std::numeric_limits<t_ilm_uint>::min(),
                                  std::numeric_limits<t_ilm_uint>::min()};
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
        t_ilm_uint surface_pos[2] = {1, 1};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surface_pos));
        ASSERT_EQ(surface_pos[0],
                  surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(surface_pos[1],
                  surfaces_allocated[i].surfaceProperties.sourceY);
    }

    // Set Orientations of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.orientation = orientation[ 3 - (i % 4) ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                            surfaces_allocated[i].surfaceProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &orientation_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.orientation,
                  orientation_rtn);
    }

    // Set position of layers to minimum
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {std::numeric_limits<t_ilm_uint>::min(),
                                   std::numeric_limits<t_ilm_uint>::min()};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check position of layers set correctly
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {1, 1};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }

    // Change things again

    // Set Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.orientation = orientation[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                                          layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &orientation_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.orientation,
                  orientation_rtn);
    }

    // Set position of surfaces to minimum
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {std::numeric_limits<t_ilm_uint>::max(),
                                  std::numeric_limits<t_ilm_uint>::min()};
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
        t_ilm_uint surface_pos[2] = {1, 1};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surface_pos));
        ASSERT_EQ(surface_pos[0],
                  surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(surface_pos[1],
                  surfaces_allocated[i].surfaceProperties.sourceY);
    }

    // Set Orientations of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.orientation = orientation[ 3 - (i % 4) ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                            surfaces_allocated[i].surfaceProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &orientation_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.orientation,
                  orientation_rtn);
    }

    // Set position of layers to minimum
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {std::numeric_limits<t_ilm_uint>::max(),
                                   std::numeric_limits<t_ilm_uint>::min()};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check position of layers set correctly
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {1, 1};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }

    // Change things again

    // Set Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.orientation = orientation[ 3 - (i % 4) ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                                          layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &orientation_rtn));
        EXPECT_EQ(layers_allocated[i].layerProperties.orientation,
                  orientation_rtn);
    }

    // Set position of surfaces to minimum
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint surf_pos[2] = {std::numeric_limits<t_ilm_uint>::min(),
                                  std::numeric_limits<t_ilm_uint>::max()};
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
        t_ilm_uint surface_pos[2] = {1, 1};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetPosition(surfaces_allocated[i].returnedSurfaceId,
                                         surface_pos));
        ASSERT_EQ(surface_pos[0],
                  surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(surface_pos[1],
                  surfaces_allocated[i].surfaceProperties.sourceY);
    }

    // Set Orientations of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        surfaces_allocated[i].surfaceProperties.orientation = orientation[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                            surfaces_allocated[i].surfaceProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check orientation of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm orientation
        e_ilmOrientation orientation_rtn;
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &orientation_rtn));
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.orientation,
                  orientation_rtn);
    }

    // Set position of layers to minimum
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {std::numeric_limits<t_ilm_uint>::min(),
                                   std::numeric_limits<t_ilm_uint>::max()};
        layers_allocated[i].layerProperties.sourceX = layer_pos[0];
        layers_allocated[i].layerProperties.sourceY = layer_pos[1];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check position of layers set correctly
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_uint layer_pos[2] = {1, 1};

        // Confirm position
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetPosition(layers_allocated[i].layerId,
                                       layer_pos));
        ASSERT_EQ(layer_pos[0], layers_allocated[i].layerProperties.sourceX);
        ASSERT_EQ(layer_pos[1], layers_allocated[i].layerProperties.sourceY);
    }

    uint num_surfaces = surfaces_allocated.size();

    // Loop through surfaces and remove
    for (uint i = 0; i < num_surfaces; i++)
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
                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetPosition(surfaceIDs[j], posreturned));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceX,
                          posreturned[0]);
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.sourceY,
                          posreturned[1]);

                // Confirm orientation
                e_ilmOrientation orientation_rtn;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetOrientation(surfaces_allocated[index].returnedSurfaceId,
                                                    &orientation_rtn));
                EXPECT_EQ(surfaces_allocated[index].surfaceProperties.orientation,
                          orientation_rtn);
            }
        }

        surfaceIDs.clear();
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
                    t_ilm_uint posreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetPosition(layerIDs[j], posreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.sourceX,
                              posreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.sourceY,
                              posreturned[1]);

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

                    // expect callback to have been called
                    assertNoCallbackIsCalled();
                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
}
