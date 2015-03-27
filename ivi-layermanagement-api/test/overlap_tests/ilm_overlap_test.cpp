/***************************************************************************
 *
 * Copyright 2015 Codethink Ltd
 * Copyright 2010-2014 BMW Car IT GmbH
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
#include <algorithm>
#include <time.h>


#include "../TestBase.h"

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

/* Tests with callbacks
 * For each test first set the global variables to point to where parameters of the callbacks are supposed to be placed.
 */
static pthread_mutex_t notificationMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  waiterVariable = PTHREAD_COND_INITIALIZER;
static int timesCalled=0;

struct surface_def {
    t_ilm_surface requestedSurfaceId;
    t_ilm_surface returnedSurfaceId;
    ilmSurfaceProperties surfaceProperties;
    bool notificationState;
};

struct layer_def {
    t_ilm_layer layerId;
    ilmLayerProperties layerProperties;
    bool notificationState;
    std::vector<t_ilm_surface> surfacesOnLayer;
};

class IlmOverlapTest : public TestBase, public ::testing::Test {

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
    std::vector<t_ilm_layer> layer_render_order;
    std::vector<t_ilm_surface> surface_render_order;
    std::vector<t_ilm_uint> v_screenID;
    static std::vector<void(IlmOverlapTest::*)(void)> vectorOfTests;
    static std::vector<std::string> vectorOfTestNames;
    static std::vector<void(IlmOverlapTest::*)(void)> vectorOfTestSet;

    static const uint no_formats = 7;
    static const uint no_surfaces = 7;
    static const uint no_layers = 4;
    static const uint no_orientations = 4;

    static const ilmPixelFormat pixelFormats[no_formats];

    static const e_ilmOrientation orientation[no_orientations];

    static bool randomize;
    static int no_iterations;
    static std::string configurationFileName;

    IlmOverlapTest() : TestBase()
    {
        if (initTests())
        {
            initSurfaces(__LINE__);
            initLayers(__LINE__);
            srand(time(NULL));
        }
    }

    ~IlmOverlapTest()
    {
        removeAll(__LINE__);
    }

    void SetUp()
    {
            vectorOfTests.push_back(&IlmOverlapTest::IlmOverlapTest_ilm_overlapGetPropertiesOfSurface);
            vectorOfTestNames.push_back("IlmOverlapTest_ilm_overlapGetPropertiesOfSurface");
            vectorOfTests.push_back(&IlmOverlapTest::IlmOverlapTest_ilm_overlapSurfaceGetDimension);
            vectorOfTestNames.push_back("IlmOverlapTest_ilm_overlapSurfaceGetDimension");
            vectorOfTests.push_back(&IlmOverlapTest::IlmOverlapTest_ilm_overlapSurfaceGetVisibility);
            vectorOfTestNames.push_back("IlmOverlapTest_ilm_overlapSurfaceGetVisibility");
    }

    void TearDown()
    {

    }

    bool initTests()
    {
        bool retval = (ILM_SUCCESS == ilm_initWithNativedisplay((t_ilm_nativedisplay)wlDisplay));

        mask = static_cast<t_ilm_notification_mask>(0);
        surface = INVALID_ID;
        timesCalled=0;

        callbackLayerId = INVALID_ID;
        callbackSurfaceId = INVALID_ID;
        return (retval);
    }

    void initSurfaces(uint lineNumber =__LINE__)
    {
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
                                         pixelFormats[i % no_formats],
                                         &(surface->returnedSurfaceId)));
            surfaces_allocated.push_back(*surface);
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }

        // Preset surface parameters
        for (uint i = 0; i < surfaces_allocated.size(); i++)
        {
            // Set dimensions of surfaces
            t_ilm_uint surf_dim[2] = {surfaces_allocated[i].surfaceProperties.origSourceWidth,
                                      surfaces_allocated[i].surfaceProperties.origSourceHeight};

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetDimension(surfaces_allocated[i].returnedSurfaceId,
                                              surf_dim));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set position of surfaces
            t_ilm_uint surf_pos[2] = {20 + (i * 5), 40 + (i * 5)};
            surfaces_allocated[i].surfaceProperties.sourceX = surf_pos[0];
            surfaces_allocated[i].surfaceProperties.sourceY = surf_pos[1];
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetPosition(surfaces_allocated[i].returnedSurfaceId,
                                             surf_pos));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set Orientations of surfaces
            surfaces_allocated[i].surfaceProperties.orientation
                = orientation[ i % no_orientations ];
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                                surfaces_allocated[i].surfaceProperties.orientation));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set Opacity of surfaces
            surfaces_allocated[i].surfaceProperties.opacity = 1.0;

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetOpacity(surfaces_allocated[i].returnedSurfaceId,
                                            surfaces_allocated[i].surfaceProperties.opacity));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set Visibility
            surfaces_allocated[i].surfaceProperties.visibility = ILM_TRUE;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetVisibility(surfaces_allocated[i].returnedSurfaceId,
                      surfaces_allocated[i].surfaceProperties.visibility));

            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set source rectangle
            surfaces_allocated[i].surfaceProperties.sourceX = 0;
            surfaces_allocated[i].surfaceProperties.sourceY = 0;
            surfaces_allocated[i].surfaceProperties.sourceWidth
                = surfaces_allocated[i].surfaceProperties.origSourceWidth;
            surfaces_allocated[i].surfaceProperties.sourceHeight
                = surfaces_allocated[i].surfaceProperties.origSourceHeight;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetSourceRectangle(surfaces_allocated[i].returnedSurfaceId,
                                                    surfaces_allocated[i].surfaceProperties.sourceX,
                                                    surfaces_allocated[i].surfaceProperties.sourceY,
                                                    surfaces_allocated[i].surfaceProperties.sourceWidth,
                                                    surfaces_allocated[i].surfaceProperties.sourceHeight));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set destination rectangle
            surfaces_allocated[i].surfaceProperties.destX
                = surfaces_allocated[i].surfaceProperties.sourceX;
            surfaces_allocated[i].surfaceProperties.destY
                = surfaces_allocated[i].surfaceProperties.sourceY;
            surfaces_allocated[i].surfaceProperties.destWidth
                = surfaces_allocated[i].surfaceProperties.sourceWidth;
            surfaces_allocated[i].surfaceProperties.destHeight
                = surfaces_allocated[i].surfaceProperties.sourceHeight;

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetDestinationRectangle(surfaces_allocated[i].returnedSurfaceId,
                                                         surfaces_allocated[i].surfaceProperties.destX,
                                                         surfaces_allocated[i].surfaceProperties.destY,
                                                         surfaces_allocated[i].surfaceProperties.destWidth,
                                                         surfaces_allocated[i].surfaceProperties.destHeight));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set notificationState
            surfaces_allocated[i].notificationState = false;
        }
    }

    void initLayers(uint lineNumber =__LINE__)
    {
        t_ilm_layer idRenderOrder[no_layers];

        // Create layers
        for (uint i = 0; i < no_layers; i++)
        {
            layer_def * layer = new layer_def;
            layer->layerId = getLayer();
            idRenderOrder[i] = layer->layerId;
            layer_render_order.push_back(layer->layerId);
            layer->layerProperties.origSourceWidth = 200 * (i + 1);
            layer->layerProperties.origSourceHeight = 240 * (i + 1);

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerCreateWithDimension(&(layer->layerId),
                                                   layer->layerProperties.origSourceWidth,
                                                   layer->layerProperties.origSourceHeight));
            layers_allocated.push_back(*layer);
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }


        for (uint i = 0; i < layers_allocated.size(); i++)
        {
            // Set position of layers check callback
            t_ilm_uint layer_pos[2] = {0, 0};
            callbackLayerId = layers_allocated[i].layerId;
            layers_allocated[i].layerProperties.sourceX = layer_pos[0];
            layers_allocated[i].layerProperties.sourceY = layer_pos[1];
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetPosition(layers_allocated[i].layerId,
                                           layer_pos));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set Orientations of layer
            layers_allocated[i].layerProperties.orientation
                = orientation[ i % no_orientations ];
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetOrientation(layers_allocated[i].layerId,
                                              layers_allocated[i].layerProperties.orientation));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set Opacity of layers
            layers_allocated[i].layerProperties.opacity = 1.0;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetOpacity(layers_allocated[i].layerId,
                      layers_allocated[i].layerProperties.opacity));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Change visibility
            layers_allocated[i].layerProperties.visibility = ILM_TRUE;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetVisibility(layers_allocated[i].layerId,
                                         layers_allocated[i].layerProperties.visibility));
            EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set source rectangle
            callbackLayerId = layers_allocated[i].layerId;

            layers_allocated[i].layerProperties.sourceX = 0;
            layers_allocated[i].layerProperties.sourceY = 0;
            layers_allocated[i].layerProperties.sourceWidth
                = layers_allocated[i].layerProperties.origSourceWidth;
            layers_allocated[i].layerProperties.sourceHeight
                = layers_allocated[i].layerProperties.origSourceHeight;

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetSourceRectangle(layers_allocated[i].layerId,
                                                  layers_allocated[i].layerProperties.sourceX,
                                                  layers_allocated[i].layerProperties.sourceY,
                                                  layers_allocated[i].layerProperties.sourceWidth,
                                                  layers_allocated[i].layerProperties.sourceHeight));

            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set destination rectangle of layers
            layers_allocated[i].layerProperties.destX = 0;
            layers_allocated[i].layerProperties.destY = 0;
            layers_allocated[i].layerProperties.destWidth
                = layers_allocated[i].layerProperties.sourceWidth;
            layers_allocated[i].layerProperties.destHeight
                = layers_allocated[i].layerProperties.sourceHeight;

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerSetDestinationRectangle(layers_allocated[i].layerId,
                                                       layers_allocated[i].layerProperties.destX,
                                                       layers_allocated[i].layerProperties.destY,
                                                       layers_allocated[i].layerProperties.destWidth,
                                                       layers_allocated[i].layerProperties.destHeight));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set notificationState
            layers_allocated[i].notificationState = false;
        }

        // Pre-set render order
        {
            t_ilm_uint numberOfScreens = 0;
            t_ilm_uint* screenIDs = NULL;
            const uint no_layers = layers_allocated.size();

            // Try to get screen IDs using null pointer for numberOfScreens
            ASSERT_EQ(ILM_FAILED, ilm_getScreenIDs(NULL, &screenIDs));

            // Try to get screen IDs using valid pointer for numberOfScreens
            ASSERT_EQ(ILM_SUCCESS, ilm_getScreenIDs(&numberOfScreens, &screenIDs));

            //Pre-set screen ID list.
            v_screenID.assign(screenIDs, screenIDs + numberOfScreens);
            free(screenIDs);

            EXPECT_TRUE(numberOfScreens>0);

            if (numberOfScreens > 0)
            {
                t_ilm_display screen = v_screenID[0];
                ilmScreenProperties screenProperties;

                ASSERT_EQ(ILM_SUCCESS,
                          ilm_displaySetRenderOrder(screen,
                                                    idRenderOrder,
                                                    layers_allocated.size()));
            }
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

                layers_allocated[i].surfacesOnLayer.push_back(surfaces_allocated[j].returnedSurfaceId);
            }
        }
    }

    void removeAll(uint lineNumber = __LINE__)
    {
        // set default values
        callbackLayerId = INVALID_ID;
        t_ilm_layer* layers = NULL;
        t_ilm_int numLayer=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&numLayer, &layers));
        for (t_ilm_int i=0; i<numLayer; i++)
        {
            if ( ( layers[i] >= getStartLayerId() )
                   && ( layers[i] <= getEndLayerId() ) )
            {
                EXPECT_EQ( ILM_SUCCESS, ilm_layerRemove( layers[i] ) );
            }
        }

        layers_allocated.clear();

        free(layers);

        t_ilm_surface* surfaces = NULL;
        t_ilm_int numSurfaces=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&numSurfaces, &surfaces));
        for (t_ilm_int i=0; i<numSurfaces; i++)
        {
            if ( ( surfaces[i] >= getStartSurfaceId() )
                    && ( layers[i] <= getEndSurfaceId() ) )
            {
                EXPECT_EQ( ILM_SUCCESS, ilm_surfaceRemove( surfaces[i] ) );
            }
        }

        free(surfaces);

        surfaces_allocated.clear();

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        ASSERT_EQ(ILM_SUCCESS, ilm_destroy());
    }

    static void assertCallbackcalled(int numberOfExpectedCalls=1, uint lineNumber =__LINE__){
        static struct timespec theTime;
        bool retval = false;
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

        EXPECT_TRUE(callbackSurfaceId == (unsigned)INVALID_ID || callbackSurfaceId == surface);
        callbackSurfaceId = surface;
        mask |= (unsigned)m;
        timesCalled++;

        pthread_cond_signal( &waiterVariable );
    }

    void IlmOverlapTest_ilm_overlapGetPropertiesOfSurface()
    {
        std::cout << "Running: " << __FUNCTION__ << std::endl;

        for (uint i = 0; i < surfaces_allocated.size(); i++)
        {
            ilmSurfaceProperties returnValue;
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                                 &returnValue))
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;

            // Check opacity
            EXPECT_NEAR(surfaces_allocated[i].surfaceProperties.opacity,
                        returnValue.opacity, 0.01)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;

            // Check source values
            ASSERT_EQ(returnValue.sourceX,
                      surfaces_allocated[i].surfaceProperties.sourceX)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
            ASSERT_EQ(returnValue.sourceY,
                      surfaces_allocated[i].surfaceProperties.sourceY)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
            ASSERT_EQ(returnValue.sourceWidth,
                      surfaces_allocated[i].surfaceProperties.sourceWidth)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
            ASSERT_EQ(returnValue.sourceHeight,
                      surfaces_allocated[i].surfaceProperties.sourceHeight)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;

            // Check destination values
            ASSERT_EQ(returnValue.sourceX,
                      surfaces_allocated[i].surfaceProperties.destX)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
            ASSERT_EQ(returnValue.sourceY,
                      surfaces_allocated[i].surfaceProperties.destY)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
            ASSERT_EQ(returnValue.sourceWidth,
                      surfaces_allocated[i].surfaceProperties.destWidth)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
            ASSERT_EQ(returnValue.sourceHeight,
                      surfaces_allocated[i].surfaceProperties.destHeight)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;

            // Check orientation value
            ASSERT_EQ(returnValue.orientation,
                      surfaces_allocated[i].surfaceProperties.orientation)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;

            // Check visibility value
            ASSERT_EQ(returnValue.visibility,
                      surfaces_allocated[i].surfaceProperties.visibility)
                              << "Surface Id: "
                              << surfaces_allocated[i].returnedSurfaceId
                              << std::endl;
        }
    }

    void IlmOverlapTest_ilm_overlapSurfaceGetDimension()
    {
        std::cout << "Running: " << __FUNCTION__ << std::endl;

        for (uint i = 0; i < surfaces_allocated.size(); i++)
        {
             t_ilm_uint dimreturned[2] = {0, 0};
             EXPECT_EQ(ILM_SUCCESS,
                       ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId,
                                               dimreturned))
                                       << "Surface Id: "
                                       << surfaces_allocated[i].returnedSurfaceId;

             EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth,
                       dimreturned[0]) << "Surface Id: "
                                       << surfaces_allocated[i].returnedSurfaceId;

             EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight,
                       dimreturned[1]) << "Surface Id: "
                                       << surfaces_allocated[i].returnedSurfaceId;
        }
    }

    void IlmOverlapTest_ilm_overlapSurfaceGetVisibility()
    {
        std::cout << "Running: " << __FUNCTION__ << std::endl;

        for (uint i = 0; i < surfaces_allocated.size(); i++)
        {
             // Confirm visibility of each surface
             t_ilm_bool visibility_rtn;
             ASSERT_EQ(ILM_SUCCESS,
                       ilm_surfaceGetVisibility(surfaces_allocated[i].returnedSurfaceId,
                                                &visibility_rtn)) << "Surface: "
                                                                  << surfaces_allocated[i].returnedSurfaceId;
             EXPECT_EQ(surfaces_allocated[i].surfaceProperties.visibility,
                       visibility_rtn)
                       << "Surface: "  << surfaces_allocated[i].returnedSurfaceId;
        }
    }
};

// Pointers where to put received values for current Test
t_ilm_layer IlmOverlapTest::callbackLayerId;
t_ilm_surface IlmOverlapTest::callbackSurfaceId;
struct ilmLayerProperties IlmOverlapTest::LayerProperties;
unsigned int IlmOverlapTest::mask;
t_ilm_surface IlmOverlapTest::surface;
ilmSurfaceProperties IlmOverlapTest::SurfaceProperties;

const ilmPixelFormat IlmOverlapTest::pixelFormats[no_formats] = {ILM_PIXELFORMAT_RGBA_4444,
                                                                 ILM_PIXELFORMAT_RGBA_5551,
                                                                 ILM_PIXELFORMAT_RGBA_6661,
                                                                 ILM_PIXELFORMAT_RGBA_8888,
                                                                 ILM_PIXELFORMAT_RGB_565,
                                                                 ILM_PIXELFORMAT_RGB_888,
                                                                 ILM_PIXELFORMAT_R_8};

const e_ilmOrientation IlmOverlapTest::orientation[no_orientations] = {ILM_ZERO,
                                                                       ILM_NINETY,
                                                                       ILM_ONEHUNDREDEIGHTY,
                                                                       ILM_TWOHUNDREDSEVENTY};

std::vector<void(IlmOverlapTest::*)(void)> IlmOverlapTest::vectorOfTests;
std::vector<std::string> IlmOverlapTest::vectorOfTestNames;
std::vector<void(IlmOverlapTest::*)(void)> IlmOverlapTest::vectorOfTestSet;
bool IlmOverlapTest::randomize;
int IlmOverlapTest::no_iterations;
std::string configurationFileName;

TEST_F(IlmOverlapTest, ilm_overlapRun)
{
    uint count = 0;

    // Note for infinity if no_iterations goes negative
    // the while loop cases will never end.

    // Check if a valid test set has been configured
    if (vectorOfTestSet.size() == 0)
    {
        // Check for random selection from all tests
        if (randomize)
        {
            // Check if the number of expected iterations have been reached
            while (count != no_iterations)
            {
                // Access a random test function
                (this->*vectorOfTests[rand() % vectorOfTests.size()])();
                count++;
            }
        }
        else
        {
            // Iterate through all tests a number of times
            while (count != no_iterations)
            {
                // Iterate through all tests
                for (uint i = 0; i < vectorOfTests.size(); i++)
                {
                    // Access next test in sequence
                    (this->*vectorOfTests[i])();
                }

                // Increment iteration counter
                count++;
            }
        }
    }
    else
    {
        // Check for random selection from specified test set
        if (randomize)
        {
            // Check if the required number of random selections is met.
            while (count != no_iterations)
            {
                // Accees a random test function from the requested set
                (this->*vectorOfTestSet[rand() % vectorOfTestSet.size()])();
                count++;
            }
        }
        else
        {
            // Iterate through all tests in the selected set a number of times
            while (count != no_iterations)
            {
                // Iterate through all tests
                for (uint i = 0; i < vectorOfTestSet.size(); i++)
                {
                    // Access next test in sequence
                    (this->*vectorOfTestSet[i])();
                }

                // Increment iteration counter
                count++;
            }
        }
    }
}
