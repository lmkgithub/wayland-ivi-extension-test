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

    static const uint no_formats = 7;
    static const uint no_surfaces = 7;
    static const uint no_layers = 4;

    IlmOverlapTest() : TestBase()
    {
        if (initTests())
        {
            initSurfaces(__LINE__);
            initLayers(__LINE__);
        }
    }

    ~IlmOverlapTest()
    {
        removeAll(__LINE__);
    }

    void SetUp()
    {

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

        ilmPixelFormat pixelFormats[no_formats] = {ILM_PIXELFORMAT_RGBA_4444,
                                               ILM_PIXELFORMAT_RGBA_5551,
                                               ILM_PIXELFORMAT_RGBA_6661,
                                               ILM_PIXELFORMAT_RGBA_8888,
                                               ILM_PIXELFORMAT_RGB_565,
                                               ILM_PIXELFORMAT_RGB_888,
                                               ILM_PIXELFORMAT_R_8};

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
            surfaces_allocated[i].surfaceProperties.orientation = orientation[ i % 4 ];
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_surfaceSetOrientation(surfaces_allocated[i].returnedSurfaceId,
                                                surfaces_allocated[i].surfaceProperties.orientation));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

            // Set Opacity of surfaces
            surfaces_allocated[i].surfaceProperties.opacity
                = 1.0;
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
            surfaces_allocated[i].surfaceProperties.sourceX
                = 0;
            surfaces_allocated[i].surfaceProperties.sourceY
                = 0;
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

        ilmPixelFormat pixelFormats[no_formats] = {ILM_PIXELFORMAT_RGBA_4444,
                                               ILM_PIXELFORMAT_RGBA_5551,
                                               ILM_PIXELFORMAT_RGBA_6661,
                                               ILM_PIXELFORMAT_RGBA_8888,
                                               ILM_PIXELFORMAT_RGB_565,
                                               ILM_PIXELFORMAT_RGB_888,
                                               ILM_PIXELFORMAT_R_8};

        e_ilmOrientation orientation[4] = {ILM_ZERO, ILM_NINETY,
                                       ILM_ONEHUNDREDEIGHTY,
                                       ILM_TWOHUNDREDSEVENTY};

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
            layers_allocated[i].layerProperties.orientation = orientation[ i % 4 ];
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
            surfaces_allocated[i].notificationState = false;
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

//        EXPECT_TRUE(callbackSurfaceId == (unsigned)INVALID_ID || callbackSurfaceId == surface);
        callbackSurfaceId = surface;
        mask |= (unsigned)m;
        timesCalled++;

        pthread_cond_signal( &waiterVariable );
    }
};

// Pointers where to put received values for current Test
t_ilm_layer IlmOverlapTest::callbackLayerId;
t_ilm_surface IlmOverlapTest::callbackSurfaceId;
struct ilmLayerProperties IlmOverlapTest::LayerProperties;
unsigned int IlmOverlapTest::mask;
t_ilm_surface IlmOverlapTest::surface;
ilmSurfaceProperties IlmOverlapTest::SurfaceProperties;

TEST_F(IlmOverlapTest, ilm_overlapGetPropertiesOfSurface)
{
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmSurfaceProperties returnValue;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId, &returnValue));

        // Check opacity
        EXPECT_NEAR(surfaces_allocated[i].surfaceProperties.opacity,
                    returnValue.opacity,
                    0.01);

        // Check source values
        ASSERT_EQ(returnValue.sourceX,
                  surfaces_allocated[i].surfaceProperties.sourceX);
        ASSERT_EQ(returnValue.sourceY,
                  surfaces_allocated[i].surfaceProperties.sourceY);
        ASSERT_EQ(returnValue.sourceWidth,
                  surfaces_allocated[i].surfaceProperties.sourceWidth);
        ASSERT_EQ(returnValue.sourceHeight,
                  surfaces_allocated[i].surfaceProperties.sourceHeight);

        // Check destination values
        ASSERT_EQ(returnValue.sourceX,
                  surfaces_allocated[i].surfaceProperties.destX);
        ASSERT_EQ(returnValue.sourceY,
                  surfaces_allocated[i].surfaceProperties.destY);
        ASSERT_EQ(returnValue.sourceWidth,
                  surfaces_allocated[i].surfaceProperties.destWidth);
        ASSERT_EQ(returnValue.sourceHeight,
                  surfaces_allocated[i].surfaceProperties.destHeight);

        // Check orientation value
        ASSERT_EQ(returnValue.orientation,
                  surfaces_allocated[i].surfaceProperties.orientation);

        // Check visibility value
        ASSERT_EQ(returnValue.visibility,
                  surfaces_allocated[i].surfaceProperties.visibility);
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceGetDimension)
{
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        t_ilm_uint dimreturned[2] = {0, 0};
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetDimension(surfaces_allocated[i].returnedSurfaceId, dimreturned));

        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceWidth,
                  dimreturned[0]);
        EXPECT_EQ(surfaces_allocated[i].surfaceProperties.origSourceHeight,
                  dimreturned[1]);
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceGetVisibility)
{
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
         // Confirm visibility of each surface
         t_ilm_bool visibility_rtn;
         ASSERT_EQ(ILM_SUCCESS,
                   ilm_surfaceGetVisibility(surfaces_allocated[i].returnedSurfaceId,
                                            &visibility_rtn));
         EXPECT_EQ(surfaces_allocated[i].surfaceProperties.visibility,
                   visibility_rtn)
                   << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                   << ", Visibility Expected: " << surfaces_allocated[i].surfaceProperties.visibility
                   << ", Visibility Got: " << visibility_rtn << std::endl;

    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceSetSourceRectangle)
{
    // Confirm source rectangles before change
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm source rectangle for each surface
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceX,
                  surfaceProperties.sourceX)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceX expected: " << surfaces_allocated[i].surfaceProperties.sourceX
                  << ", sourceX got: " << surfaceProperties.sourceX << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceY,
                  surfaceProperties.sourceY)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceY expected: " << surfaces_allocated[i].surfaceProperties.sourceY
                  << ", sourceY got: " << surfaceProperties.sourceY << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceWidth,
                  surfaceProperties.sourceWidth)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceWidth expected: " << surfaces_allocated[i].surfaceProperties.sourceWidth
                  << ", sourceWidth got: " << surfaceProperties.sourceWidth << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceHeight,
                  surfaceProperties.sourceHeight)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceHeight expected: " << surfaces_allocated[i].surfaceProperties.sourceHeight
                  << ", sourceHeight got: " << surfaceProperties.sourceHeight << std::endl;
    }

    if (surfaces_allocated.size() > 0)
    {
        // Set random surface index
        uint random_surface = rand() % surfaces_allocated.size();

        // Create random values
        t_ilm_uint random_sourceX = rand();
        t_ilm_uint random_sourceY = rand();
        t_ilm_uint random_sourceWidth = rand();
        t_ilm_uint random_sourceHeight = rand();

        // Set callback
        callbackSurfaceId = surfaces_allocated[random_surface].returnedSurfaceId;

        // Set source rectangle
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetSourceRectangle(surfaces_allocated[random_surface].returnedSurfaceId,
                                                random_sourceX,
                                                random_sourceY,
                                                random_sourceWidth,
                                                random_sourceHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Update stored orientation for surface
        surfaces_allocated[random_surface].surfaceProperties.sourceX = random_sourceX;
        surfaces_allocated[random_surface].surfaceProperties.sourceY = random_sourceY;
        surfaces_allocated[random_surface].surfaceProperties.sourceWidth = random_sourceWidth; 
        surfaces_allocated[random_surface].surfaceProperties.sourceHeight = random_sourceHeight;

        // Check notification state if set
        if (surfaces_allocated[random_surface].notificationState)
        {
            assertCallbackcalled();
        }
    }

    // Confirm all source rectangles after change
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm source rectangle for each surface
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceX,
                  surfaceProperties.sourceX)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceX expected: " << surfaces_allocated[i].surfaceProperties.sourceX
                  << ", sourceX got: " << surfaceProperties.sourceX << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceY,
                  surfaceProperties.sourceY)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceY expected: " << surfaces_allocated[i].surfaceProperties.sourceY
                  << ", sourceY got: " << surfaceProperties.sourceY << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceWidth,
                  surfaceProperties.sourceWidth)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceWidth expected: " << surfaces_allocated[i].surfaceProperties.sourceWidth
                  << ", sourceWidth got: " << surfaceProperties.sourceWidth << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceHeight,
                  surfaceProperties.sourceHeight)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceHeight expected: " << surfaces_allocated[i].surfaceProperties.sourceHeight
                  << ", sourceHeight got: " << surfaceProperties.sourceHeight << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapGetPropertiesOfLayer)
{
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmLayerProperties returnValue;
        // Check Layer source properties
        EXPECT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers_allocated[i].layerId, &returnValue));

        // Check opacity
        EXPECT_NEAR(layers_allocated[i].layerProperties.opacity,
                    returnValue.opacity,
                    0.01);

        // Check source values
        EXPECT_EQ(returnValue.sourceX,
                  layers_allocated[i].layerProperties.sourceX);
        EXPECT_EQ(returnValue.sourceY,
                  layers_allocated[i].layerProperties.sourceY);
        EXPECT_EQ(returnValue.sourceWidth,
                  layers_allocated[i].layerProperties.sourceWidth);
        EXPECT_EQ(returnValue.sourceHeight,
                  layers_allocated[i].layerProperties.sourceHeight);

        // Check destination values
        EXPECT_EQ(returnValue.sourceX,
                  layers_allocated[i].layerProperties.destX);
        EXPECT_EQ(returnValue.sourceY,
                  layers_allocated[i].layerProperties.destY);
        EXPECT_EQ(returnValue.sourceWidth,
                  layers_allocated[i].layerProperties.destWidth);
        EXPECT_EQ(returnValue.sourceHeight,
                  layers_allocated[i].layerProperties.destHeight);

        // Check orientation value
        EXPECT_EQ(returnValue.orientation,
                  layers_allocated[i].layerProperties.orientation);

        // Check visibility value
        EXPECT_EQ(returnValue.visibility,
                  layers_allocated[i].layerProperties.visibility);
    }
}

TEST_F(IlmOverlapTest, ilm_overlapGetScreenIDs)
{
    t_ilm_uint numberOfScreens = 0;
    t_ilm_uint* screenIDs;

    // Try to get screen IDs using valid pointer for numberOfScreens
    ASSERT_EQ(ILM_SUCCESS, ilm_getScreenIDs(&numberOfScreens, &screenIDs));

    v_screenID.clear();
    v_screenID.assign(screenIDs, screenIDs + numberOfScreens);
    free(screenIDs);

    EXPECT_TRUE(numberOfScreens > 0);
}

TEST_F(IlmOverlapTest, ilm_overlapGetLayerIDs)
{
    t_ilm_int length;
    t_ilm_layer* IDs;
    std::vector<t_ilm_layer> layerIDs;

    // Get layers
    ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));
    layerIDs.assign(IDs, IDs + length);
    free(IDs);

    // Loop through expected layers and confirm they are present
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        EXPECT_NE(std::find(layerIDs.begin(),
                            layerIDs.end(),
                            layers_allocated[i].layerId)
                            , layerIDs.end());
    }
}

TEST_F(IlmOverlapTest, ilm_overlapGetLayerIDsOnScreen)
{
    t_ilm_layer* idGotRenderOrder;
    t_ilm_int    length = 0;
    std::vector<t_ilm_layer> layerIDs;

    // This assumes only 1 screen for the moment
    ASSERT_EQ(ILM_SUCCESS,
              ilm_getLayerIDsOnScreen(v_screenID[0],
                                      &length,
                                      &idGotRenderOrder));

    layerIDs.assign(idGotRenderOrder, idGotRenderOrder + length);
    free(idGotRenderOrder);

    // Loop through expected layers and confirm they are present
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ASSERT_NE(std::find(layerIDs.begin(),
                            layerIDs.end(),
                            layers_allocated[i].layerId)
                            , layerIDs.end());
    }

    // Check that the number of layers got match those expected.
    if (layerIDs.size() != layer_render_order.size())
    {
        ASSERT_EQ(true, false) << "No of layers retrieved don't match expected"
                               << "Got: " << layerIDs.size() << ", "
                               << "Expected: " << layer_render_order.size()
                               << std::endl;
    }

    // Check got render order matches that expected
    for (uint i = 0; i < layerIDs.size(); i++)
    {
        ASSERT_EQ(layer_render_order[i], layerIDs[i]);
    }

    // Clean-up
    layerIDs.clear();
}

TEST_F(IlmOverlapTest, ilm_overlapGetSurfaceIDs)
{
    t_ilm_int length;
    t_ilm_surface* IDs;
    std::vector<t_ilm_surface> surfaceIDs;

    // Get surfaces
    ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));
    surfaceIDs.assign(IDs, IDs + length);
    free(IDs);

    // Loop through expected surfaces and confirm they are present
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        EXPECT_NE(std::find(surfaceIDs.begin(),
                            surfaceIDs.end(),
                            surfaces_allocated[i].returnedSurfaceId)
                            , surfaceIDs.end());
    }
}

TEST_F(IlmOverlapTest, ilm_overlapGetSurfaceIDsOnLayer)
{
    t_ilm_int length;
    t_ilm_uint* IDs;
    std::vector<t_ilm_surface> surfaceIDs;

    // Loop through expected layers and confirm they are present
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getSurfaceIDsOnLayer(layers_allocated[i].layerId, &length, &IDs));
        surfaceIDs.assign(IDs, IDs + length);
        free(IDs);

        // Check to make sure that all the expected surfaces are assigned.
        if (surfaceIDs.size() == layers_allocated[i].surfacesOnLayer.size())
        {
            for (uint j = 0; j < surfaceIDs.size(); j++)
            {
                 EXPECT_NE(std::find(layers_allocated[i].surfacesOnLayer.begin(),
                           layers_allocated[i].surfacesOnLayer.end(),
                           surfaceIDs[j]),
                           layers_allocated[i].surfacesOnLayer.end());
            }
        }
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceGetOrientation)
{
    // Check Orientations of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &returned));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.orientation, returned)
            << " - Surface: " << surfaces_allocated[i].returnedSurfaceId
            << ", Orientation Expected: " << surfaces_allocated[i].surfaceProperties.orientation
            << ", Orientation Got: " << returned << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceSetOrientation)
{
    e_ilmOrientation orientation[4] = {ILM_ZERO, ILM_NINETY,
                                       ILM_ONEHUNDREDEIGHTY,
                                       ILM_TWOHUNDREDSEVENTY};

    // Check Orientations of surfaces
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &returned));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.orientation, returned)
            << " - Surface: " << surfaces_allocated[i].returnedSurfaceId
            << ", Orientation Expected: " << surfaces_allocated[i].surfaceProperties.orientation
            << ", Orientation Got: " << returned << std::endl;
    }

    // Pick random surface and change
    uint random_surface = rand() % surfaces_allocated.size();
    e_ilmOrientation random_orientation = orientation[rand() % 4];
    callbackSurfaceId = surfaces_allocated[random_surface].returnedSurfaceId;

    if (surfaces_allocated.size() > 0)
    {
        // Set a random surface a random orientation
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetOrientation(surfaces_allocated[random_surface].returnedSurfaceId,
                                            random_orientation));

        // Update stored orientation for surface
        surfaces_allocated[random_surface].surfaceProperties.orientation = random_orientation;

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Check notification state if set
        if (surfaces_allocated[random_surface].notificationState)
        {
            assertCallbackcalled();
        }
    }

    // Check Orientations of surfaces again
    // Make sure changing one hasn't modified others
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceGetOrientation(surfaces_allocated[i].returnedSurfaceId,
                  &returned));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.orientation, returned)
            << " - Surface: " << surfaces_allocated[i].returnedSurfaceId
            << ", Orientation Expected: " << surfaces_allocated[i].surfaceProperties.orientation
            << ", Orientation Got: " << returned << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceSetVisibility)
{
    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
         // Confirm visibility of each surface
         t_ilm_bool visibility_rtn;
         ASSERT_EQ(ILM_SUCCESS,
                   ilm_surfaceGetVisibility(surfaces_allocated[i].returnedSurfaceId,
                                            &visibility_rtn));
         EXPECT_EQ(surfaces_allocated[i].surfaceProperties.visibility,
                   visibility_rtn)
                   << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                   << ", Visibility Expected: " << surfaces_allocated[i].surfaceProperties.visibility
                   << ", Visibility Got: " << visibility_rtn << std::endl;

    }

    if (surfaces_allocated.size() > 0)
    {
        // Set random surface index
        uint random_surface = rand() % surfaces_allocated.size();

        // Create random value for visibility
        t_ilm_bool random_visibility = visibility[rand() % 2];

        // Set callback
        callbackSurfaceId = surfaces_allocated[random_surface].returnedSurfaceId;

        // Set visibility
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetVisibility(surfaces_allocated[random_surface].returnedSurfaceId,
                                           random_visibility));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Update stored orientation for surface
        surfaces_allocated[random_surface].surfaceProperties.visibility = random_visibility;

        // Check notification state if set
        if (surfaces_allocated[random_surface].notificationState)
        {
            assertCallbackcalled();
        }
    }

    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
         // Confirm visibility of each surface after change
         t_ilm_bool visibility_rtn;
         ASSERT_EQ(ILM_SUCCESS,
                   ilm_surfaceGetVisibility(surfaces_allocated[i].returnedSurfaceId,
                                            &visibility_rtn));
         EXPECT_EQ(surfaces_allocated[i].surfaceProperties.visibility,
                   visibility_rtn)
                   << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                   << ", Visibility Expected: " << surfaces_allocated[i].surfaceProperties.visibility
                   << ", Visibility Got: " << visibility_rtn << std::endl;

    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceGetDestinationRectangle)
{
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm source rectangle
        ilmSurfaceProperties surfaceProperties;
        // Confirm destination rectangle
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destX,
                  surfaceProperties.destX)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destX expected: " << surfaces_allocated[i].surfaceProperties.destX
                  << ", destX got: " << surfaceProperties.destX << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destY,
                  surfaceProperties.destY)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destY expected: " << surfaces_allocated[i].surfaceProperties.destY
                  << ", destY got: " << surfaceProperties.destY << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destWidth,
                  surfaceProperties.destWidth)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destWidth expected: " << surfaces_allocated[i].surfaceProperties.destWidth
                  << ", destWidth got: " << surfaceProperties.destWidth << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destHeight,
                  surfaceProperties.destHeight)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destHeight expected: " << surfaces_allocated[i].surfaceProperties.destHeight
                  << ", destHeight got: " << surfaceProperties.destHeight << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceSetDestinationRectangle)
{
    // Confirm destination rectangles before change
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm destination rectangle for each surface
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destX,
                  surfaceProperties.destX)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destX expected: " << surfaces_allocated[i].surfaceProperties.destX
                  << ", destX got: " << surfaceProperties.destX << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destY,
                  surfaceProperties.destY)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destY expected: " << surfaces_allocated[i].surfaceProperties.destY
                  << ", destY got: " << surfaceProperties.destY << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destWidth,
                  surfaceProperties.destWidth)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destWidth expected: " << surfaces_allocated[i].surfaceProperties.destWidth
                  << ", destWidth got: " << surfaceProperties.destWidth << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destHeight,
                  surfaceProperties.destHeight)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destHeight expected: " << surfaces_allocated[i].surfaceProperties.destHeight
                  << ", destHeight got: " << surfaceProperties.destHeight << std::endl;
    }

    if (surfaces_allocated.size() > 0)
    {
        // Set random surface index
        uint random_surface = rand() % surfaces_allocated.size();

        // Create random values
        t_ilm_uint random_destX = rand();
        t_ilm_uint random_destY = rand();
        t_ilm_uint random_destWidth = rand();
        t_ilm_uint random_destHeight = rand();

        // Set callback
        callbackSurfaceId = surfaces_allocated[random_surface].returnedSurfaceId;

        // Set destination rectangle
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetDestinationRectangle(surfaces_allocated[random_surface].returnedSurfaceId,
                                                     random_destX,
                                                     random_destY,
                                                     random_destWidth,
                                                     random_destHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Update stored orientation for surface
        surfaces_allocated[random_surface].surfaceProperties.destX = random_destX;
        surfaces_allocated[random_surface].surfaceProperties.destY = random_destY;
        surfaces_allocated[random_surface].surfaceProperties.destWidth = random_destWidth; 
        surfaces_allocated[random_surface].surfaceProperties.destHeight = random_destHeight;

        // Check notification state if set
        if (surfaces_allocated[random_surface].notificationState)
        {
            assertCallbackcalled();
        }
    }

    // Confirm all destination rectangles after change
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm destination rectangle for each surface
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destX,
                  surfaceProperties.destX)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destX expected: " << surfaces_allocated[i].surfaceProperties.destX
                  << ", destX got: " << surfaceProperties.destX << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destY,
                  surfaceProperties.destY)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destY expected: " << surfaces_allocated[i].surfaceProperties.destY
                  << ", destY got: " << surfaceProperties.destY << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destWidth,
                  surfaceProperties.destWidth)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destWidth expected: " << surfaces_allocated[i].surfaceProperties.destWidth
                  << ", destWidth got: " << surfaceProperties.destWidth << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destHeight,
                  surfaceProperties.destHeight)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", destHeight expected: " << surfaces_allocated[i].surfaceProperties.destHeight
                  << ", destHeight got: " << surfaceProperties.destHeight << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapSurfaceGetSourceRectangle)
{
    for (uint i = 0; i < surfaces_allocated.size(); i++)
    {
        // Confirm source rectangle
        ilmSurfaceProperties surfaceProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_allocated[i].returnedSurfaceId,
                                             &surfaceProperties));
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceX,
                  surfaceProperties.sourceX)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceX expected: " << surfaces_allocated[i].surfaceProperties.sourceX
                  << ", sourceX got: " << surfaceProperties.sourceX << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceY,
                  surfaceProperties.sourceY)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceY expected: " << surfaces_allocated[i].surfaceProperties.sourceY
                  << ", sourceY got: " << surfaceProperties.sourceY << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.destWidth,
                  surfaceProperties.destWidth)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceWidth expected: " << surfaces_allocated[i].surfaceProperties.sourceWidth
                  << ", sourceWidth got: " << surfaceProperties.sourceWidth << std::endl;
        ASSERT_EQ(surfaces_allocated[i].surfaceProperties.sourceHeight,
                  surfaceProperties.sourceHeight)
                  << "Surface: "  << surfaces_allocated[i].returnedSurfaceId
                  << ", sourceHeight expected: " << surfaces_allocated[i].surfaceProperties.sourceHeight
                  << ", sourceHeight got: " << surfaceProperties.sourceHeight << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapLayerGetOrientation)
{
    // Check Orientations of surfaces
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &returned));
        ASSERT_EQ(layers_allocated[i].layerProperties.orientation, returned)
            << " - Layer: " << layers_allocated[i].layerId
            << ", Orientation Expected: " << layers_allocated[i].layerProperties.orientation
            << ", Orientation Got: " << returned << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapLayerSetOrientation)
{
    e_ilmOrientation orientation[4] = {ILM_ZERO, ILM_NINETY,
                                       ILM_ONEHUNDREDEIGHTY,
                                       ILM_TWOHUNDREDSEVENTY};

    // Check Orientations of layers
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &returned));
        ASSERT_EQ(layers_allocated[i].layerProperties.orientation, returned)
            << " - Layer: " << layers_allocated[i].layerId
            << ", Orientation Expected: " << layers_allocated[i].layerProperties.orientation
            << ", Orientation Got: " << returned << std::endl;
    }

    // Pick random layer and change
    uint random_layer = rand() % layers_allocated.size();
    e_ilmOrientation random_orientation = orientation[rand() % 4];
    callbackLayerId = layers_allocated[random_layer].layerId;

    if (layers_allocated.size() > 0)
    {
        // Set a random layer a random orientation
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetOrientation(layers_allocated[random_layer].layerId,
                                          random_orientation));

        // Update stored orientation for layer
        layers_allocated[random_layer].layerProperties.orientation
            = random_orientation;

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Check notification state if set
        if (layers_allocated[random_layer].notificationState)
        {
            assertCallbackcalled();
        }
    }

    // Check Orientations of layers again
    // Make sure changing one hasn't modified others
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerGetOrientation(layers_allocated[i].layerId,
                  &returned));
        ASSERT_EQ(layers_allocated[i].layerProperties.orientation, returned)
            << " - Layers: " << layers_allocated[i].layerId
            << ", Orientation Expected: " << layers_allocated[i].layerProperties.orientation
            << ", Orientation Got: " << returned << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapLayerGetVisibility)
{
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
         // Confirm visibility of each surface
         t_ilm_bool visibility_rtn;
         ASSERT_EQ(ILM_SUCCESS,
                   ilm_layerGetVisibility(layers_allocated[i].layerId,
                                          &visibility_rtn));
         EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                   visibility_rtn)
                   << "Layer: "  << layers_allocated[i].layerId
                   << ", Visibility Expected: " << layers_allocated[i].layerProperties.visibility
                   << ", Visibility Got: " << visibility_rtn << std::endl;

    }
}

TEST_F(IlmOverlapTest, ilm_overlapLayerSetVisibility)
{
    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    for (uint i = 0; i < layers_allocated.size(); i++)
    {
         // Confirm visibility of each layer
         t_ilm_bool visibility_rtn;
         ASSERT_EQ(ILM_SUCCESS,
                   ilm_layerGetVisibility(layers_allocated[i].layerId,
                                            &visibility_rtn));
         EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                   visibility_rtn)
                   << "Layer: "  << layers_allocated[i].layerId
                   << ", Visibility Expected: " << layers_allocated[i].layerProperties.visibility
                   << ", Visibility Got: " << visibility_rtn << std::endl;

    }

    if (layers_allocated.size() > 0)
    {
        // Set random layer index
        uint random_layer = rand() % layers_allocated.size();

        // Create random value for visibility
        t_ilm_bool random_visibility = visibility[rand() % 2];

        // Set callback
        callbackLayerId = layers_allocated[random_layer].layerId;

        // Set visibility
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetVisibility(layers_allocated[random_layer].layerId,
                                         random_visibility));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Update stored orientation for layer
        layers_allocated[random_layer].layerProperties.visibility = random_visibility;

        // Check notification state if set
        if (layers_allocated[random_layer].notificationState)
        {
            assertCallbackcalled();
        }
    }

    for (uint i = 0; i < layers_allocated.size(); i++)
    {
         // Confirm visibility of each layer after change
         t_ilm_bool visibility_rtn;
         ASSERT_EQ(ILM_SUCCESS,
                   ilm_layerGetVisibility(layers_allocated[i].layerId,
                                            &visibility_rtn));
         EXPECT_EQ(layers_allocated[i].layerProperties.visibility,
                   visibility_rtn)
                   << "Layer: "  << layers_allocated[i].layerId
                   << ", Visibility Expected: " << layers_allocated[i].layerProperties.visibility
                   << ", Visibility Got: " << visibility_rtn << std::endl;

    }
}

TEST_F(IlmOverlapTest, ilm_overlapLayerGetDestinationRectangle)
{
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm source rectangle
        ilmLayerProperties layerProperties;
        // Confirm destination rectangle
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers_allocated[i].layerId,
                                           &layerProperties));
        ASSERT_EQ(layers_allocated[i].layerProperties.destX,
                  layerProperties.destX)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destX expected: " << layers_allocated[i].layerProperties.destX
                  << ", destX got: " << layerProperties.destX << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destY,
                  layerProperties.destY)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destY expected: " << layers_allocated[i].layerProperties.destY
                  << ", destY got: " << layerProperties.destY << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destWidth,
                  layerProperties.destWidth)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destWidth expected: " << layers_allocated[i].layerProperties.destWidth
                  << ", destWidth got: " << layerProperties.destWidth << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destHeight,
                  layerProperties.destHeight)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destHeight expected: " << layers_allocated[i].layerProperties.destHeight
                  << ", destHeight got: " << layerProperties.destHeight << std::endl;
    }
}

TEST_F(IlmOverlapTest, ilm_overlapLayerSetDestinationRectangle)
{
    // Confirm destination rectangles before change
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm destination rectangle for each layer
        ilmLayerProperties layerProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers_allocated[i].layerId,
                                             &layerProperties));
        ASSERT_EQ(layers_allocated[i].layerProperties.destX,
                  layerProperties.destX)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destX expected: " << layers_allocated[i].layerProperties.destX
                  << ", destX got: " << layerProperties.destX << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destY,
                  layerProperties.destY)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destY expected: " << layers_allocated[i].layerProperties.destY
                  << ", destY got: " << layerProperties.destY << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destWidth,
                  layerProperties.destWidth)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destWidth expected: " << layers_allocated[i].layerProperties.destWidth
                  << ", destWidth got: " << layerProperties.destWidth << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destHeight,
                  layerProperties.destHeight)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destHeight expected: " << layers_allocated[i].layerProperties.destHeight
                  << ", destHeight got: " << layerProperties.destHeight << std::endl;
    }

    if (layers_allocated.size() > 0)
    {
        // Set random layer index
        uint random_layer = rand() % layers_allocated.size();

        // Create random values
        t_ilm_uint random_destX = rand();
        t_ilm_uint random_destY = rand();
        t_ilm_uint random_destWidth = rand();
        t_ilm_uint random_destHeight = rand();

        // Set callback
        callbackLayerId = layers_allocated[random_layer].layerId;

        // Set destination rectangle
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetDestinationRectangle(layers_allocated[random_layer].layerId,
                                                     random_destX,
                                                     random_destY,
                                                     random_destWidth,
                                                     random_destHeight));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Update stored orientation for layer
        layers_allocated[random_layer].layerProperties.destX = random_destX;
        layers_allocated[random_layer].layerProperties.destY = random_destY;
        layers_allocated[random_layer].layerProperties.destWidth = random_destWidth; 
        layers_allocated[random_layer].layerProperties.destHeight = random_destHeight;

        // Check notification state if set
        if (layers_allocated[random_layer].notificationState)
        {
            assertCallbackcalled();
        }
    }

    // Confirm all destination rectangles after change
    for (uint i = 0; i < layers_allocated.size(); i++)
    {
        // Confirm destination rectangle for each layer
        ilmLayerProperties layerProperties;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers_allocated[i].layerId,
                                             &layerProperties));
        ASSERT_EQ(layers_allocated[i].layerProperties.destX,
                  layerProperties.destX)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destX expected: " << layers_allocated[i].layerProperties.destX
                  << ", destX got: " << layerProperties.destX << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destY,
                  layerProperties.destY)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destY expected: " << layers_allocated[i].layerProperties.destY
                  << ", destY got: " << layerProperties.destY << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destWidth,
                  layerProperties.destWidth)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destWidth expected: " << layers_allocated[i].layerProperties.destWidth
                  << ", destWidth got: " << layerProperties.destWidth << std::endl;
        ASSERT_EQ(layers_allocated[i].layerProperties.destHeight,
                  layerProperties.destHeight)
                  << "Layer: "  << layers_allocated[i].layerId
                  << ", destHeight expected: " << layers_allocated[i].layerProperties.destHeight
                  << ", destHeight got: " << layerProperties.destHeight << std::endl;
    }
}
