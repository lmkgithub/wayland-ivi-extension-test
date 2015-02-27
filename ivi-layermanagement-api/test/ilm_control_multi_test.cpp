/***************************************************************************
 *
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

#include "TestBase.h"

extern "C" {
    #include "ilm_client.h"
    #include "ilm_control.h"
}

template <typename T>
bool contains(T const *actual, size_t as, T expected)
{
   for (unsigned i = 0; i < as; i++)
      if (actual[i] == expected)
         return true;
   return false;
}

class IlmCommandMultiTest : public TestBase, public ::testing::Test {
public:
    void SetUp()
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_initWithNativedisplay((t_ilm_nativedisplay)wlDisplay));
    }

    void TearDown()
    {
        t_ilm_layer* layers = NULL;
        t_ilm_int numLayers=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&numLayers, &layers));
        for (t_ilm_int i=0; i<numLayers; i++)
        {
            EXPECT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        };
        free(layers);

        t_ilm_surface* surfaces = NULL;
        t_ilm_int numSurfaces=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&numSurfaces, &surfaces));
        for (t_ilm_int i=0; i<numSurfaces; i++)
        {
            EXPECT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[i]));
        };
        free(surfaces);

        EXPECT_EQ(ILM_SUCCESS, ilm_commitChanges());
        EXPECT_EQ(ILM_SUCCESS, ilm_destroy());
    }
};

TEST_F(IlmCommandMultiTest, multi_SetGetSurfaceDimension) {
    t_ilm_surface surface_start;
    t_ilm_layer layer_start;
    uint no_surfaces = 4;
    uint no_layers = 2;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_dim[no_layers][2];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 15 * (i + 1);
        surf_dim[i][1] = 25 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS, 
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i], 
                                     surf_dim[i][0], surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();
        layer_dim[i][0] = 200 * (i + 1);
        layer_dim[i][1] = 240 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0],
                                               layer_dim[i][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[i], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add 1st set of surfaces to 1st layer
    for (int i = 0; i < (no_surfaces / no_layers); i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_layerAddSurface(layers[0], surfaces[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add 2nd set of surfaces to 2nd layer
    for (int i = (no_surfaces / no_layers); i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_layerAddSurface(layers[1], surfaces[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for surfaces and compare against those originally set
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces[i], dimreturned));
        EXPECT_EQ(surf_dim[i][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[i][1], dimreturned[1]); 
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int i = 0; i < length; i++)
        {
            uint index = no_surfaces;
            for (int j = 0; j < no_surfaces; j++)
            {
                if (IDs[i] == surfaces_dest[j]) index = j;
            }

            t_ilm_uint dimreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[i], dimreturned));
            EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
            EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
        }

        free(IDs);
    }


    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layer
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check dimensions 
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint dimreturned[2] = {0, 0};
            uint layer_index = IDs[j] - layers[0];
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned)); 
            EXPECT_EQ(layer_dim[layer_index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[layer_index][1], dimreturned[1]);
        }

        free(IDs);
    }
}


TEST_F(IlmCommandMultiTest, multi_SetGetLayerDimension) {
    uint no_surfaces = 16;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_dim[no_layers][2];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 * (i + 1);
        surf_dim[i][1] =  15 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0], surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();
        layer_dim[i][0] = 100 * (i + 1);
        layer_dim[i][1] = 120 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0],
                                               layer_dim[i][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for layers and compare against those originally set
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layers[i], dimreturned));
        EXPECT_EQ(layer_dim[i][0], dimreturned[0]);
        EXPECT_EQ(layer_dim[i][1], dimreturned[1]);
    }

    // Change dimensions of layers
    for (uint i = 0; i < no_layers; i++)
    {
        surf_dim[i][0] = 12 * (i + 1);
        surf_dim[i][1] =  17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetDimension(layers[i], layer_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for layers and compare against modified values
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(layers[i], dimreturned));
        EXPECT_EQ(layer_dim[i][0], dimreturned[0]);
        EXPECT_EQ(layer_dim[i][1], dimreturned[1]);
    }

    // Set dimensions of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[i], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }
    
    // Add surfaces to layers
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = (i * (no_surfaces / no_layers));
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers[i],
                                          surfaces[j]));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }

    }

    // Get dimensions for surfaces and compare against those originally set
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces[i], dimreturned));
        EXPECT_EQ(surf_dim[i][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[i][1], dimreturned[1]);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));


        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }
        free(IDs);
    }


    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check dimensions 
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint dimreturned[2] = {0, 0};
            uint layer_index = IDs[j] - layers[0];
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[layer_index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[layer_index][1], dimreturned[1]);
        }

        free(IDs);
    }

}

TEST_F(IlmCommandMultiTest, multi_SetGetSurfacePosition) {

    t_ilm_uint surf_pos_start[2] = {15, 25};
    t_ilm_uint layer_pos_start[2] = {25, 40};
    t_ilm_surface surface_start = 44;
    t_ilm_layer layer_start = 4306;
    uint no_surfaces = 4;
    uint no_layers = 2;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_dim[no_layers][2];
    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];


    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0] = 10 + (i * 5),
                                     surf_dim[i][1] = 10 + (i * 5),
                                     ILM_PIXELFORMAT_RGBA_8888, &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0] = 200 * (i + 1),
                                               layer_dim[i][1] = 240 * (i + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = surf_pos_start[0] + (i * 5);
        surf_pos[i][1] = surf_pos_start[1] + (i * 5);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = layer_pos_start[0] + (i * 10);
        layer_pos[i][1] = layer_pos_start[1] + (i * 10);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }


    ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    // Check and clear any previous surfaces on layer
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[i], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = (i * (no_surfaces / no_layers));
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers[i],
                                          surfaces[j]));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }

    }

    // Get positions for surfaces and compare against those originally set
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint posreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaces[i], posreturned));
        EXPECT_EQ(surf_pos[i][0], posreturned[0]);
        EXPECT_EQ(surf_pos[i][1], posreturned[1]);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm positions are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layer
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check position
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));
            for (int k = 0; k < no_layers; k++)
            {
                if (IDs[j] == layers[k])
                {
                    EXPECT_EQ(layer_pos[k][0], posreturned[0]);
                    EXPECT_EQ(layer_pos[k][1], posreturned[1]);
                }
            }
        }

        free(IDs);
    }

}

TEST_F(IlmCommandMultiTest, multi_SetGetLayerPosition) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 20 + (i * 15);
        surf_dim[i][1] = 20 + (i * 15);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0] = 105 * (i + 1),
                                               layer_dim[i][1] = 125 * (i + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (uint i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 5 * (i * 20);
        surf_pos[i][1] = 7 * (i * 25);

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check and clear any previous surfaces
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[i], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers
    for (int i = 0; i < no_layers; i++)
    {
        for (int j = (i * (no_surfaces / no_layers));
             j < ((i + 1) * (no_surfaces / no_layers));
             j++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers[i],
                                          surfaces[j]));
            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
        }

    }

    // Get dimensions for surfaces and compare against those originally set
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces_dest[i], dimreturned));
        EXPECT_EQ(surf_dim[i][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[i][1], dimreturned[1]);
    }

    // Get positions for surfaces and compare against those originally set
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint posreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaces[i], posreturned));
        EXPECT_EQ(surf_pos[i][0], posreturned[0]);
        EXPECT_EQ(surf_pos[i][1], posreturned[1]);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions & positions
        // are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            uint layer_index = IDs[j] - layers[0];
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[layer_index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[layer_index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[layer_index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[layer_index][1], dimreturned[1]);

        }

        free(IDs);
    }

}


TEST_F(IlmCommandMultiTest, multi_SetGetSurfaceOrientation) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];
    e_ilmOrientation orientation[4] = {ILM_ZERO, ILM_NINETY,
                                       ILM_ONEHUNDREDEIGHTY,
                                       ILM_TWOHUNDREDSEVENTY};

    e_ilmOrientation surface_orientations[no_surfaces];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 + (i * 35);
        surf_dim[i][1] = 15 + (i * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for surfaces and compare against those originally set
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces_dest[i], dimreturned));
        EXPECT_EQ(surf_dim[i][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[i][1], dimreturned[1]);
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();
        layer_dim[i][0] = 85 * (i + 1);
        layer_dim[i][1] = 75 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0],
                                               layer_dim[i][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
    }

    // Set Orientations of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_orientations[i] = orientation[ i % 4 ];
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetOrientation(surfaces_dest[i], surface_orientations[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Orientations of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetOrientation(surfaces_dest[i], &returned));
        ASSERT_EQ(surface_orientations[i], returned);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                e_ilmOrientation orientation_returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetOrientation(IDs[j],
                                                    &orientation_returned));
                ASSERT_EQ(surface_orientations[index], orientation_returned);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            uint layer_index = IDs[j] - layers[0];
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[layer_index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[layer_index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[layer_index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[layer_index][1], dimreturned[1]);

        }

        free(IDs);
    }
}

TEST_F(IlmCommandMultiTest, multi_SetGetLayerOrientation) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];
    e_ilmOrientation orientation[4] = {ILM_ZERO, ILM_NINETY,
                                       ILM_ONEHUNDREDEIGHTY,
                                       ILM_TWOHUNDREDSEVENTY};

    e_ilmOrientation layer_orientations[no_surfaces];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 + (i * 35);
        surf_dim[i][1] = 15 + (i * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();
        layer_dim[i][0] = 85 * (i + 1);
        layer_dim[i][1] = 75 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0],
                                               layer_dim[i][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Orientations of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_orientations[i] = orientation[ i % 4];
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetOrientation(layers[i], layer_orientations[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Orientations of layers
    for (int i = 0; i < no_layers; i++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOrientation(layers[i], &returned));
        ASSERT_EQ(layer_orientations[i], returned);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                e_ilmOrientation orientation_returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            e_ilmOrientation orientation_returned;
            uint layer_index = IDs[j] - layers[0];

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[layer_index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[layer_index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[layer_index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[layer_index][1], dimreturned[1]);

            EXPECT_EQ(ILM_SUCCESS,
                      ilm_layerGetOrientation(IDs[j],
                                                &orientation_returned));
            ASSERT_EQ(layer_orientations[layer_index], orientation_returned);
        }

        free(IDs);
    }
}

TEST_F(IlmCommandMultiTest, multi_SetGetLayerSurfaceOpacity) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    t_ilm_float surface_opacitys[no_surfaces];
    t_ilm_float layer_opacitys[no_layers];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 + (i * 35);
        surf_dim[i][1] = 15 + (i * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();
        layer_dim[i][0] = 85 * (i + 1);
        layer_dim[i][1] = 75 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0],
                                               layer_dim[i][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Opacity of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_opacitys[i] = 0.1 + (i * 0.2);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetOpacity(surfaces_dest[i], surface_opacitys[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Opacity of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_opacitys[i] = 0.05 + (i * 0.15);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetOpacity(layers[i], layer_opacitys[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Opacity of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetOpacity(surfaces_dest[i], &returned));
        EXPECT_NEAR(surface_opacitys[i], returned, 0.01);
    }

    // Check Opacity of layers
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layers[i], &returned));
        EXPECT_NEAR(layer_opacitys[i], returned, 0.01);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetOpacity(surfaces_dest[index], &returned));
                EXPECT_NEAR(surface_opacitys[index], returned, 0.01);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_float returned;
            uint index = IDs[j] - layers[0];

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));
            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));
            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);

            ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(IDs[j], &returned));
            EXPECT_NEAR(layer_opacitys[index], returned, 0.01);
        }

        free(IDs);
    }
}

TEST_F(IlmCommandMultiTest, multi_SetGetSurfaceVisibility) {
    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    t_ilm_bool surface_visibility[no_surfaces];
    t_ilm_bool layer_visibility[no_layers];

    t_ilm_bool visibility[2] = {ILM_TRUE, ILM_FALSE};

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 + (i * 35);
        surf_dim[i][1] = 15 + (i * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0] = 85 * (i + 1),
                                               layer_dim[i][1] = 75 * (i + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surface_visibility[i] = visibility[ i % 2];
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetVisibility(surfaces_dest[i], surface_visibility[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_visibility[i] = visibility[ i % 2];
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetVisibility(layers[i], layer_visibility[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Visibility of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_bool returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetVisibility(surfaces_dest[i], &returned));
        ASSERT_EQ(surface_visibility[i], returned);
    }

    // Check Visibility of layers
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_bool returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetVisibility(layers[i], &returned));
        ASSERT_EQ(layer_visibility[i], returned);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                t_ilm_bool returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetVisibility(surfaces_dest[index], &returned));
                ASSERT_EQ(surface_visibility[index], returned);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_bool returned;
            t_ilm_int index = IDs[j] - layers[0];

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerGetVisibility(IDs[j], &returned));
            ASSERT_EQ(layer_visibility[index], returned);
        }

        free(IDs);
    }
}


TEST_F(IlmCommandMultiTest, multi_SetSurfaceSourceRectangle) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    ilmSurfaceProperties surfaces_Properties[no_surfaces];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 + (i * 35);
        surf_dim[i][1] = 15 + (i * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set source rectangle of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces_Properties[i].sourceX = 89 + (i * 10);
        surfaces_Properties[i].sourceY = 6538 + (i * 3);
        surfaces_Properties[i].sourceWidth = 638 + (i * 7);
        surfaces_Properties[i].sourceHeight = 4 + (i * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetSourceRectangle(surfaces_dest[i],
                                                surfaces_Properties[i].sourceX,
                                                surfaces_Properties[i].sourceY,
                                                surfaces_Properties[i].sourceWidth,
                                                surfaces_Properties[i].sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm source rectangles of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ilmSurfaceProperties returnValue;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_dest[i], &returnValue));

        ASSERT_EQ(returnValue.sourceX, surfaces_Properties[i].sourceX);
        ASSERT_EQ(returnValue.sourceY, surfaces_Properties[i].sourceY);
        ASSERT_EQ(returnValue.sourceWidth, surfaces_Properties[i].sourceWidth);
        ASSERT_EQ(returnValue.sourceHeight, surfaces_Properties[i].sourceHeight);
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();
        layer_dim[i][0] = 85 * (i + 1);
        layer_dim[i][0] = 75 * (i + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0],
                                               layer_dim[i][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfSurface(IDs[j], &returnValue));
                ASSERT_EQ(returnValue.sourceX, surfaces_Properties[index].sourceX);
                ASSERT_EQ(returnValue.sourceY, surfaces_Properties[index].sourceY);
                ASSERT_EQ(returnValue.sourceWidth, surfaces_Properties[index].sourceWidth);
                ASSERT_EQ(returnValue.sourceHeight, surfaces_Properties[index].sourceHeight);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[j] - layers[0];

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
        }

        free(IDs);
    }
}

TEST_F(IlmCommandMultiTest, multi_SetLayerSourceRectangle) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    ilmLayerProperties layers_Properties[no_layers];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 10 + (i * 35);
        surf_dim[i][1] = 15 + (i * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0] = 85 * (i + 1),
                                               layer_dim[i][1] = 75 * (i + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 12 * (i + 1);
        layer_pos[i][1] = 17 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set source rectangle of layers
    for (int i = 0; i < no_layers; i++)
    {
        layers_Properties[i].sourceX = 79 + (i * 10);
        layers_Properties[i].sourceY = 6438 + (i * 3);
        layers_Properties[i].sourceWidth = 538 + (i * 7);
        layers_Properties[i].sourceHeight = 2 + (i * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers[i],
                                              layers_Properties[i].sourceX,
                                              layers_Properties[i].sourceY,
                                              layers_Properties[i].sourceWidth,
                                              layers_Properties[i].sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm source rectangles of layers
    for (int i = 0; i < no_layers; i++)
    {
        ilmLayerProperties returnValue;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers[i], &returnValue));

        ASSERT_EQ(returnValue.sourceX, layers_Properties[i].sourceX);
        ASSERT_EQ(returnValue.sourceY, layers_Properties[i].sourceY);
        ASSERT_EQ(returnValue.sourceWidth, layers_Properties[i].sourceWidth);
        ASSERT_EQ(returnValue.sourceHeight, layers_Properties[i].sourceHeight);
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[j] - layers[0];
            ilmLayerProperties returnValue;

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfLayer(IDs[j], &returnValue));
            ASSERT_EQ(returnValue.sourceX, layers_Properties[index].sourceX);
            ASSERT_EQ(returnValue.sourceY, layers_Properties[index].sourceY);
            ASSERT_EQ(returnValue.sourceWidth, layers_Properties[index].sourceWidth);
            ASSERT_EQ(returnValue.sourceHeight, layers_Properties[index].sourceHeight);

        }

        free(IDs);
    }
}

TEST_F(IlmCommandMultiTest, ilm_multiTakeSurfaceLayerScreenshots) {

    uint no_surfaces = 4;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    const char* outputFile = "/tmp/test.bmp";

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 8 + (i * 15);
        surf_dim[i][1] = 9 + (i * 35);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0] = 25 * (i + 1),
                                               layer_dim[i][1] = 85 * (i + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 14 * (i + 1);
        layer_pos[i][1] = 16 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 12 * (i * 32);
        surf_pos[i][1] = 13 * (i * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop round and take layer shots
    for (int i = 0; i < no_layers; i++)
    {
        // make sure the file is not there before
        FILE* f = fopen(outputFile, "r");
        if (f!=NULL){
            fclose(f);
            int result = remove(outputFile);
            ASSERT_EQ(0, result);
        }

        ASSERT_EQ(ILM_SUCCESS, ilm_takeLayerScreenshot(outputFile, layers[i]));

        sleep(1);
        f = fopen(outputFile, "r");
        ASSERT_TRUE(f!=NULL);
        fclose(f);
        remove(outputFile);
    }

    // Loop round and take surface shots
    for (int i = 0; i < no_surfaces; i++)
    {
        // make sure the file is not there before
        FILE* f = fopen(outputFile, "r");
        if (f!=NULL){
            fclose(f);
            int result = remove(outputFile);
            ASSERT_EQ(0, result);
        }

        ASSERT_EQ(ILM_SUCCESS, ilm_takeSurfaceScreenshot(outputFile, surfaces_dest[i]));

        sleep(1);
        f = fopen(outputFile, "r");
        ASSERT_TRUE(f!=NULL);
        fclose(f);
        remove(outputFile);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[j] - layers[0];
            ilmLayerProperties returnValue;

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfLayer(IDs[j], &returnValue));

        }

        free(IDs);
    }
}


TEST_F(IlmCommandMultiTest, ilm_multiLayerSurfaceGetPixelformat) {
    t_ilm_uint surface1=0;
    t_ilm_uint surface2=1;
    t_ilm_uint surface3=2;
    t_ilm_uint surface4=3;
    t_ilm_uint surface5=4;
    t_ilm_uint surface6=5;
    t_ilm_uint surface7=6;

    t_ilm_surface surface_start = 68;
    t_ilm_layer layer_start = 4330;
    uint no_surfaces = 8;
    uint no_layers = 4;
    t_ilm_surface surfaces[no_surfaces];
    t_ilm_surface surfaces_dest[no_surfaces];
    t_ilm_layer layers[no_layers];

    t_ilm_uint surf_pos[no_surfaces][2];
    t_ilm_uint surf_dim[no_surfaces][2];
    t_ilm_uint layer_pos[no_layers][2];
    t_ilm_uint layer_dim[no_layers][2];

    e_ilmPixelFormat formats[8] = {ILM_PIXELFORMAT_R_8, ILM_PIXELFORMAT_RGB_888,
                                   ILM_PIXELFORMAT_RGBA_8888,
                                   ILM_PIXELFORMAT_RGB_565,
                                   ILM_PIXELFORMAT_RGBA_5551,
                                   ILM_PIXELFORMAT_RGBA_6661,
                                   ILM_PIXELFORMAT_RGBA_4444,
                                   ILM_PIXEL_FORMAT_UNKNOWN};

    ilmPixelFormat surface_formats[no_surfaces];

    // Create surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surfaces[i] = getSurface();
        surfaces_dest[i] = surfaces[i];
        surf_dim[i][0] = 18 + (i * 12);
        surf_dim[i][1] = 29 + (i * 37);
        surface_formats[i] = formats[i % 8];

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[i],
                                     surf_dim[i][0],
                                     surf_dim[i][1],
                                     surface_formats[i],
                                     &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[i], surf_dim[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int i = 0; i < no_layers; i++)
    {
        layers[i] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[i]),
                                               layer_dim[i][0] = 32 * (i + 1),
                                               layer_dim[i][1] = 72 * (i + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int i = 0; i < no_layers; i++)
    {
        layer_pos[i][0] = 14 * (i + 1);
        layer_pos[i][1] = 16 * (i + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[i], layer_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        surf_pos[i][0] = 15 * (i * 38);
        surf_pos[i][1] = 17 * (i * 42);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[i], surf_pos[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check pixel format of surfaces
    for (int i = 0; i < no_surfaces; i++)
    {
        ilmPixelFormat pf;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetPixelformat(surfaces_dest[i], &pf));
        ASSERT_EQ(surface_formats[i], pf);
    }

    // Loop through surfaces and remove
    for (int i = 0; i < no_surfaces; i++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int j = 0; j < length; j++)
        {
            // Find the surface reference from the ID returned
            uint index = no_surfaces;
            for (int k = 0; k < no_surfaces; k++)
            {
                if (IDs[j] == surfaces_dest[k]) index = k;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                ilmPixelFormat pf_returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[j], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[j], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetPixelformat(IDs[j], &pf_returned));
                ASSERT_EQ(surface_formats[index], pf_returned);
            }
            else
            {
                ASSERT_NE(index, no_surfaces);
            }
        }

        free(IDs);
    }

    // remove the layers
    for (int i = 0; i < no_layers; i++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[i]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int j = 0; j < length; j++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[j] - layers[0];
            ilmLayerProperties returnValue;

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[j], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[j], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfLayer(IDs[j], &returnValue));

        }

        free(IDs);
    }
}
