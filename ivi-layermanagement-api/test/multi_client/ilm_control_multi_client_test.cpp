/***************************************************************************
 *
 * Copyright 2010-2014 BMW Car IT GmbH
 * Copyright (C) 2012 DENSO CORPORATION and Robert Bosch Car Multimedia Gmbh
 * Copyright (C) 2015 Codethink Ltd
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

static std::vector<t_ilm_layer> layer_callbacks;
static std::vector<t_ilm_surface> surface_callbacks;

struct surface_def {
    t_ilm_uint requestedSurfaceId;
    t_ilm_uint returnedSurfaceId;
    ilmSurfaceProperties surfaceProperties;
};

struct layer_def {
    t_ilm_uint layerId;
    ilmLayerProperties layerProperties;
};

class IlmCommandMultiClientTest : public TestBase, public ::testing::Test {

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

    IlmCommandMultiClientTest() : TestBase() {}
    ~IlmCommandMultiClientTest() {}

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
        for (t_ilm_int i=0; i < numLayer; i++)
        {
            if ( ( layers[i] >= getStartLayerId() )
                   && ( layers[i] <= getEndLayerId() ) )
            {
                EXPECT_EQ( ILM_SUCCESS, ilm_layerRemove( layers[i] ) );
            }
        };

        layers_allocated.clear();

        free(layers);

        t_ilm_surface* surfaces = NULL;
        t_ilm_int numSurfaces=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&numSurfaces, &surfaces));
        for (t_ilm_int i=0; i <numSurfaces; i++)
        {
            if ( ( surfaces[i] >= getStartSurfaceId() )
                   && ( surfaces[i] <= getEndSurfaceId() ) )
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

        layer_callbacks.push_back(layer);

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

        callbackSurfaceId = surface;
        mask |= (unsigned)m;
        timesCalled++;

        surface_callbacks.push_back(surface);

        pthread_cond_signal( &waiterVariable );
    }
};

// Pointers where to put received values for current Test
t_ilm_layer IlmCommandMultiClientTest::callbackLayerId;
t_ilm_surface IlmCommandMultiClientTest::callbackSurfaceId;
struct ilmLayerProperties IlmCommandMultiClientTest::LayerProperties;
unsigned int IlmCommandMultiClientTest::mask;
t_ilm_surface IlmCommandMultiClientTest::surface;
ilmSurfaceProperties IlmCommandMultiClientTest::SurfaceProperties;

/*****************************************************************************
 * This test allows simultaneous multiple client testing. Each test creates  *
 * a number of surfaces (with notifications) and layers and then sets a      *
 * number of parameters (dimension, opacity, source rectangle) checks that   *
 * notifications have been called when required. The surfaces are then added *
 * to layers. The parameters previously set are checked to confirm that      *
 * adding surfaces to layers does not change underlying parameter values.    *
 * Surfaces and layers are removed checking on each occasion that the values *
 * of other surfaces and layers are unchanged by the process of removal.     *
 *                                                                           *
 * Each test runner can be configured with separate surface and layer sets   *
 * (to prevent clashes during multi-client testing) through the command line *
 * interface. An example multi-client script is provided to illustrate this. *
 *****************************************************************************/
TEST_F(IlmCommandMultiClientTest, multi_Test) {
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

    // Add notifications to surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // add notification
        ilmErrorTypes status = ilm_surfaceAddNotification(surfaces_allocated[i].returnedSurfaceId,&SurfaceCallbackFunction);
        ASSERT_EQ(ILM_SUCCESS, status);
    }

    // Create layers
    for (uint i = 0; i < no_layers; i++)
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

    // Set dimensions of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        callbackSurfaceId = surfaces_allocated[i].returnedSurfaceId;
        t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                  surfaces_allocated[i].surfaceProperties.origSourceHeight};

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                          surf_dim));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Expect correct number of callback to have been called
    assertCallbackcalled(surfaces_allocated.size());

    // Expect the callbacks to tally up with the order of changes
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        EXPECT_EQ(surfaces_allocated[i].returnedSurfaceId,
                  surface_callbacks[i]);
    }

    // Clear callbacks registered
    surface_callbacks.clear();

    // Set Opacity of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        callbackSurfaceId = surfaces_allocated[i].returnedSurfaceId;
        surfaces_allocated[i].surfaceProperties.opacity = 0.05 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOpacity(surfaces_allocated[i].returnedSurfaceId,
                  surfaces_allocated[i].surfaceProperties.opacity));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Expect correct number of callback to have been called
    assertCallbackcalled(surfaces_allocated.size());

    // Expect the callbacks to tally up with the order of changes
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        EXPECT_EQ(surfaces_allocated[i].returnedSurfaceId,
                  surface_callbacks[i]);
    }

    // Clear callbacks registered
    surface_callbacks.clear();

    // Set source rectangle of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        callbackSurfaceId = surfaces_allocated[i].returnedSurfaceId;
        surfaces_allocated[i].surfaceProperties.sourceX = 79 + (i * 10);
        surfaces_allocated[i].surfaceProperties.sourceY = 6508 + (i * 3);
        surfaces_allocated[i].surfaceProperties.sourceWidth = 618 + (i * 7);
        surfaces_allocated[i].surfaceProperties.sourceHeight = 6 + (i * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetSourceRectangle(surfaces_allocated[i].returnedSurfaceId,
                                                surfaces_allocated[i].surfaceProperties.sourceX,
                                                surfaces_allocated[i].surfaceProperties.sourceY,
                                                surfaces_allocated[i].surfaceProperties.sourceWidth,
                                                surfaces_allocated[i].surfaceProperties.sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Expect correct number of callback to have been called
    assertCallbackcalled(surfaces_allocated.size());

    // Expect the callbacks to tally up with the order of changes
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        EXPECT_EQ(surfaces_allocated[i].returnedSurfaceId,
                  surface_callbacks[i]);
    }

    // Clear callbacks registered
    surface_callbacks.clear();

    // Confirm source rectangles of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmSurfaceProperties returnValue;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId, &returnValue));

        ASSERT_EQ(returnValue.sourceX, surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(returnValue.sourceY, surfaces_allocated[i].surfaceProperties.sourceY);
        ASSERT_EQ(returnValue.sourceWidth, surfaces_allocated[i].surfaceProperties.sourceWidth);
        ASSERT_EQ(returnValue.sourceHeight, surfaces_allocated[i].surfaceProperties.sourceHeight);
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

    // Check Opacity of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetOpacity(layers_allocated[i].layerId, &returned));
        EXPECT_NEAR(layers_allocated[i].layerProperties.opacity, returned, 0.01);
    }

    // Set Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        layers_allocated[i].layerProperties.orientation = orientation[i % 4];
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[i].layerId,
                                          layers_allocated[i].layerProperties.orientation));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
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
        for (int j = i * (surfaces_allocated.size() / layers_allocated.size());
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

    // Check Opacity of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOpacity(surfaces_allocated[i].returnedSurfaceId,
                                        &returned));
        EXPECT_NEAR(surfaces_allocated[i].surfaceProperties.opacity,
                    returned,
                    0.01);
    }

    // Check Opacity of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layers_allocated[i].layerId,
                                                   &returned));
        EXPECT_NEAR(layers_allocated[i].layerProperties.opacity, returned, 0.01);
    }

    // Check Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &returned));
        ASSERT_EQ(layers_allocated[i].layerProperties.orientation, returned);
    }

    // Loop through surfaces and remove
    uint num_surfaces = surfaces_allocated.size();

    std::vector<t_ilm_surface> surfaceIDCallbacksCalled;

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

                // Check Opacity of surfaces
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetOpacity(surfaceIDs[j], &returned));
                EXPECT_NEAR(surfaces_allocated[index].surfaceProperties.opacity,
                            returned, 0.01);

                // Check Surface sour properties
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfSurface(surfaceIDs[j], &returnValue));
                EXPECT_EQ(returnValue.sourceX,
                          surfaces_allocated[index].surfaceProperties.sourceX);
                EXPECT_EQ(returnValue.sourceY,
                          surfaces_allocated[index].surfaceProperties.sourceY);
                EXPECT_EQ(returnValue.sourceWidth,
                          surfaces_allocated[index].surfaceProperties.sourceWidth);
                EXPECT_EQ(returnValue.sourceHeight,
                          surfaces_allocated[index].surfaceProperties.sourceHeight);

                // Confirm that notifications are still in place on other surfaces
                // Change opacity value check then change it back. check callbacks
                surfaces_allocated[index].surfaceProperties.opacity =
                    surfaces_allocated[index].surfaceProperties.opacity - 1.0;
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceSetOpacity(surfaceIDs[j],
                          surfaces_allocated[index].surfaceProperties.opacity));
                EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());
                callbackSurfaceId = surfaceIDs[j];
                surfaceIDCallbacksCalled.push_back(surfaceIDs[j]);

                EXPECT_EQ(callbackSurfaceId,
                          surfaces_allocated[index].returnedSurfaceId);

                surfaces_allocated[index].surfaceProperties.opacity =
                    surfaces_allocated[index].surfaceProperties.opacity + 1.0;

                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceSetOpacity(surfaceIDs[j],
                          surfaces_allocated[index].surfaceProperties.opacity));
                EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());
                callbackSurfaceId = surfaceIDs[j];
                surfaceIDCallbacksCalled.push_back(surfaceIDs[j]);
                EXPECT_EQ(callbackSurfaceId,
                          surfaces_allocated[index].returnedSurfaceId);
            }
        }

        surfaceIDs.clear();
    }

    assertCallbackcalled(surfaceIDCallbacksCalled.size());

    // Check that the callbacks made match those expected.
    bool valid = (surfaceIDCallbacksCalled == surface_callbacks);

    EXPECT_EQ(valid, true);

    // Clear requested callback vector
    surfaceIDCallbacksCalled.clear();

    // Clear registered callback vector
    surface_callbacks.clear();

    surfaces_allocated.clear();

    // remove the layers
    uint num_layers = layers_allocated.size();

    for (uint i = 0; i < num_layers; i++)
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
        for (uint j = 0; j < length; j++)
        {
            uint index = num_layers;

            for (uint k = 0; k < layers_allocated.size(); k++)
            {
                if (layerIDs[j] == layers_allocated[k].layerId)
                {
                    index = k;
                    break;
                }
            }

            if (index != num_layers)
            {
                // Iterate round remaining layers and check dimensions
                for (uint s = 0; s < length; s++)
                {
                    t_ilm_uint dimreturned[2] = {0, 0};
                    e_ilmOrientation orientation_returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetDimension(layerIDs[j], dimreturned));

                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceWidth,
                              dimreturned[0]);
                    EXPECT_EQ(layers_allocated[index].layerProperties.origSourceHeight,
                              dimreturned[1]);

                    t_ilm_float returned;
                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetOpacity(layers_allocated[index].layerId,
                              &returned));
                    EXPECT_NEAR(layers_allocated[index].layerProperties.opacity,
                                returned, 0.01);

                    EXPECT_EQ(ILM_SUCCESS,
                              ilm_layerGetOrientation(layerIDs[j],
                                                      &orientation_returned));
                    EXPECT_EQ(layers_allocated[index].layerProperties.orientation,
                              orientation_returned);

                }
            }
        }

        layerIDs.clear();
    }

    layers_allocated.clear();
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
                        std::cout << "Starting Layer ID set to: "
                                  << ID << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid layer start parameter"
                                           << std::endl;
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
                        std::cout << "Starting Surface ID set to: "
                                  << ID << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid surface start parameter"
                                           << std::endl;
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
                        std::cout << "Maximum number of contiguous layers "
                                  << "set to: " << maxNumberLayers
                                  << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid maximum layer parameter"
                                           << std::endl;
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
                        std::cout << "Maximum number of contiguous surfaces "
                                  << "set to: " << maxNumberSurfaces
                                  << std::endl;
                    }
                    else
                    {
                        EXPECT_TRUE(false) << "Invalid maximum surface "
                                           << "parameter" << std::endl;
                    }
                }
            break;

            default:
            {
                 EXPECT_TRUE(false) << "Unknown parameter specified: "
                                    << argv[1][3] << std::endl;
                 EXPECT_TRUE(false) << "Options: --l<starting layer id>"
                                       << std::endl
                                    << "         --s<starting surface id>"
                                       << std::endl
                                    << "         --n<number of layers to be used>"
                                       << std::endl
                                    << "         --t<number of surfaces to be used>"
                                       << std::endl;
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
