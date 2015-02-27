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
        //print_lmc_get_scene();
        t_ilm_layer* layers = NULL;
        t_ilm_int numLayer=0;
        EXPECT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&numLayer, &layers));
        for (t_ilm_int i=0; i<numLayer; i++)
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
                                     ILM_PIXELFORMAT_RGBA_8888, &(surfaces_dest[i])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Create layers
    for (int j = 0; j < no_layers; j++)
    {
        layers[j] = getLayer();
        layer_dim[j][0] = 200 * (j + 1);
        layer_dim[j][1] = 240 * (j + 1);
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[j]),
                                               layer_dim[j][0],
                                               layer_dim[j][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint k = 0; k < no_surfaces; k++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[k], surf_dim[k]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    // Check and clear any previous surfaces
    for (int j = 0; j < no_layers; j++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[j], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add 1st set of surfaces to 1st layer
    for (int l = 0; l < no_surfaces/no_layers; l++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_layerAddSurface(layers[0], surfaces[l]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Add 2nd set of surfaces to 2nd layer
    for (int m = no_surfaces/no_layers; m < no_surfaces; m++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_layerAddSurface(layers[1], surfaces[m]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for surfaces and compare against those originally set
    for (int n = 0; n < no_surfaces; n++)
    {
        t_ilm_uint dimreturned[no_surfaces];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces[n], dimreturned));
        EXPECT_EQ(surf_dim[n][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[n][1], dimreturned[1]); 
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            t_ilm_uint dimreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
            EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
            EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
        }

        free(IDs);
    }


    // remove the layers
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layer
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check dimensions 
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint dimreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned)); 

            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][1], dimreturned[1]);

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
    for (int j = 0; j < no_layers; j++)
    {
        layers[j] = getLayer();
        layer_dim[j][0] = 100 * (j + 1);
        layer_dim[j][1] = 120 * (j + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[j]),
                                               layer_dim[j][0],
                                               layer_dim[j][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set dimensions of surfaces
    for (uint k = 0; k < no_surfaces; k++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[k], surf_dim[k]));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    // Check and clear any previous surfaces
    for (int j = 0; j < no_layers; j++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[j], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }
    
    // Add surfaces to layers
    for (int layer_count = 0; layer_count < no_layers; layer_count++)
    {
        for (int surface_count = (layer_count * (no_surfaces/no_layers));
             surface_count < ((layer_count + 1) * (no_surfaces/no_layers));
             surface_count++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers[layer_count],
                                          surfaces[surface_count]));

            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        }

    }

    // Get dimensions for surfaces and compare against those originally set
    for (int n = 0; n < no_surfaces; n++)
    {
        t_ilm_uint dimreturned[no_surfaces];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces[n], dimreturned));
        EXPECT_EQ(surf_dim[n][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[n][1], dimreturned[1]);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));


        // Loop through remaining surfaces and confirm dimensions are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check dimensions 
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint dimreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][1], dimreturned[1]);

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
    for (int j = 0; j < no_layers; j++)
    {
        layers[j] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[j]),
                                               layer_dim[j][0] = 200 * (j + 1),
                                               layer_dim[j][1] = 240 * (j + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (uint k = 0; k < no_surfaces; k++)
    {
        surf_pos[k][0] = surf_pos_start[0] + (k * 5);
        surf_pos[k][1] = surf_pos_start[1] + (k * 5);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces[k], surf_pos[k]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint l = 0; l < no_layers; l++)
    {
        layer_pos[l][0] = layer_pos_start[0] + (l * 10);
        layer_pos[l][1] = layer_pos_start[1] + (l * 10);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[l], layer_pos[l]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }


    ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    // Check and clear any previous surfaces on layer
    for (int j = 0; j < no_layers; j++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[j], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers
    for (int layer_count = 0; layer_count < no_layers; layer_count++)
    {
        for (int surface_count = (layer_count * (no_surfaces/no_layers));
             surface_count < ((layer_count + 1) * (no_surfaces/no_layers));
             surface_count++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers[layer_count],
                                          surfaces[surface_count]));

            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        }

    }

    // Get positions for surfaces and compare against those originally set
    for (int n = 0; n < no_surfaces; n++)
    {
        t_ilm_uint posreturned[no_surfaces];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaces[n], posreturned));
        EXPECT_EQ(surf_pos[n][0], posreturned[0]);
        EXPECT_EQ(surf_pos[n][1], posreturned[1]);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm positions are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layer
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check position
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));
            for (int t = 0; t < no_layers; t++)
            {
                if (IDs[s] == layers[t])
                {
                    EXPECT_EQ(layer_pos[t][0], posreturned[0]);
                    EXPECT_EQ(layer_pos[t][1], posreturned[1]);
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
    for (int j = 0; j < no_layers; j++)
    {
        layers[j] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[j]),
                                               layer_dim[j][0] = 105 * (j + 1),
                                               layer_dim[j][1] = 125 * (j + 1)));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (uint k = 0; k < no_layers; k++)
    {
        layer_pos[k][0] = 12 * (k + 1);
        layer_pos[k][1] = 17 * (k + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[k], layer_pos[k]));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (uint l = 0; l < no_surfaces; l++)
    {
        surf_pos[l][0] = 5 * (l * 20);
        surf_pos[l][1] = 7 * (l * 25);

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces[l], surf_pos[l]));

        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    // Check and clear any previous surfaces
    for (int j = 0; j < no_layers; j++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDsOnLayer(layers[j], &length, &IDs));
        free(IDs);
        ASSERT_EQ(length, 0);
    }

    // Add surfaces to layers
    for (int layer_count = 0; layer_count < no_layers; layer_count++)
    {
        for (int surface_count = (layer_count * (no_surfaces/no_layers));
             surface_count < ((layer_count + 1) * (no_surfaces/no_layers));
             surface_count++)
        {
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerAddSurface(layers[layer_count],
                                          surfaces[surface_count]));

            ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        }

    }

    // Get dimensions for surfaces and compare against those originally set
    for (int n = 0; n < no_surfaces; n++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces_dest[n], dimreturned));
        EXPECT_EQ(surf_dim[n][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[n][1], dimreturned[1]);
    }

    // Get positions for surfaces and compare against those originally set
    for (int o = 0; o < no_surfaces; o++)
    {
        t_ilm_uint posreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(surfaces[o], posreturned));
        EXPECT_EQ(surf_pos[o][0], posreturned[0]);
        EXPECT_EQ(surf_pos[o][1], posreturned[1]);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions & positions
        // are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2] = {0, 0};
                t_ilm_uint posreturned[2] = {0, 0};
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][0], posreturned[0]);
            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][1], dimreturned[1]);

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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 10 + (b * 35);
        surf_dim[b][1] = 15 + (b * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Get dimensions for surfaces and compare against those originally set
    for (int n = 0; n < no_surfaces; n++)
    {
        t_ilm_uint dimreturned[2];
        EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(surfaces_dest[n], dimreturned));
        EXPECT_EQ(surf_dim[n][0], dimreturned[0]);
        EXPECT_EQ(surf_dim[n][1], dimreturned[1]);
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();
        layer_dim[c][0] = 85 * (c + 1);
        layer_dim[c][1] = 75 * (c + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0],
                                               layer_dim[c][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 12 * (f + 1);
        layer_pos[f][1] = 17 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
    }

    // Set Orientations of surfaces
    for (int d = 0; d < no_surfaces; d++)
    {
        surface_orientations[d] = orientation[d/4];
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetOrientation(surfaces_dest[d], surface_orientations[d]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Orientations of surfaces
    for (int e = 0; e < no_surfaces; e++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetOrientation(surfaces_dest[e], &returned));
        ASSERT_EQ(surface_orientations[e], returned);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                e_ilmOrientation orientation_returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                EXPECT_EQ(ILM_SUCCESS,
                          ilm_surfaceGetOrientation(IDs[q],
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][0], posreturned[0]);
            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][1], dimreturned[1]);

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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 10 + (b * 35);
        surf_dim[b][1] = 15 + (b * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();
        layer_dim[c][0] = 85 * (c + 1);
        layer_dim[c][1] = 75 * (c + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0],
                                               layer_dim[c][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 12 * (f + 1);
        layer_pos[f][1] = 17 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Orientations of layers
    for (int d = 0; d < no_layers; d++)
    {
        layer_orientations[d] = orientation[d/4];
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetOrientation(layers[d], layer_orientations[d]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Orientations of layers
    for (int e = 0; e < no_layers; e++)
    {
        ilmOrientation returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOrientation(layers[e], &returned));
        ASSERT_EQ(layer_orientations[e], returned);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                e_ilmOrientation orientation_returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            e_ilmOrientation orientation_returned;
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][0], posreturned[0]);
            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][1], dimreturned[1]);

            EXPECT_EQ(ILM_SUCCESS,
                      ilm_layerGetOrientation(IDs[s],
                                                &orientation_returned));
            ASSERT_EQ(layer_orientations[s], orientation_returned);
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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 10 + (b * 35);
        surf_dim[b][1] = 15 + (b * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();
        layer_dim[c][0] = 85 * (c + 1);
        layer_dim[c][1] = 75 * (c + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0],
                                               layer_dim[c][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 12 * (f + 1);
        layer_pos[f][1] = 17 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Opacity of surfaces
    for (int d = 0; d < no_surfaces; d++)
    {
        surface_opacitys[d] = 0.1 + (d * 0.2);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetOpacity(surfaces_dest[d], surface_opacitys[d]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Opacity of layers
    for (int d = 0; d < no_layers; d++)
    {
        layer_opacitys[d] = 0.05 + (d * 0.15);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetOpacity(layers[d], layer_opacitys[d]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Opacity of surfaces
    for (int e = 0; e < no_surfaces; e++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetOpacity(surfaces_dest[e], &returned));
        EXPECT_NEAR(surface_opacitys[e], returned, 0.01);
    }

    // Check Opacity of layers
    for (int e = 0; e < no_surfaces; e++)
    {
        t_ilm_float returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(layers[e], &returned));
        EXPECT_NEAR(layer_opacitys[e], returned, 0.01);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                t_ilm_float returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_float returned;
            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][0], posreturned[0]);
            EXPECT_EQ(layer_pos[IDs[s] - layers[0]][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[IDs[s] - layers[0]][1], dimreturned[1]);

            ASSERT_EQ(ILM_SUCCESS, ilm_layerGetOpacity(IDs[s], &returned));
            EXPECT_NEAR(layer_opacitys[IDs[s] - layers[0]], returned, 0.01);
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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 10 + (b * 35);
        surf_dim[b][1] = 15 + (b * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0] = 85 * (c + 1),
                                               layer_dim[c][1] = 75 * (c + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 12 * (f + 1);
        layer_pos[f][1] = 17 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility of surfaces
    for (int d = 0; d < no_surfaces; d++)
    {
        surface_visibility[d] = visibility[d/2];
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetVisibility(surfaces_dest[d], surface_visibility[d]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set Visibility of layers
    for (int d = 0; d < no_layers; d++)
    {
        layer_visibility[d] = visibility[d/2];
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetVisibility(layers[d], layer_visibility[d]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check Visibility of surfaces
    for (int e = 0; e < no_surfaces; e++)
    {
        t_ilm_bool returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetVisibility(surfaces_dest[e], &returned));
        ASSERT_EQ(surface_visibility[e], returned);
    }

    // Check Visibility of layers
    for (int e = 0; e < no_surfaces; e++)
    {
        t_ilm_bool returned;
        ASSERT_EQ(ILM_SUCCESS, ilm_layerGetVisibility(layers[e], &returned));
        ASSERT_EQ(layer_visibility[e], returned);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                t_ilm_bool returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_bool returned;
            t_ilm_int index = IDs[s] - layers[0];

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);

            ASSERT_EQ(ILM_SUCCESS,
                      ilm_layerGetVisibility(IDs[s], &returned));
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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 10 + (b * 35);
        surf_dim[b][1] = 15 + (b * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set source rectangle of surfaces
    for (int d = 0; d < no_surfaces; d++)
    {
        surfaces_Properties[d].sourceX = 89 + (d * 10);
        surfaces_Properties[d].sourceY = 6538 + (d * 3);
        surfaces_Properties[d].sourceWidth = 638 + (d * 7);
        surfaces_Properties[d].sourceHeight = 4 + (d * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceSetSourceRectangle(surfaces_dest[d],
                                                surfaces_Properties[d].sourceX,
                                                surfaces_Properties[d].sourceY,
                                                surfaces_Properties[d].sourceWidth,
                                                surfaces_Properties[d].sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm source rectangles of surfaces
    for (int d = 0; d < no_surfaces; d++)
    {
        ilmSurfaceProperties returnValue;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfSurface(surfaces_dest[d], &returnValue));

        ASSERT_EQ(returnValue.sourceX, surfaces_Properties[d].sourceX);
        ASSERT_EQ(returnValue.sourceY, surfaces_Properties[d].sourceY);
        ASSERT_EQ(returnValue.sourceWidth, surfaces_Properties[d].sourceWidth);
        ASSERT_EQ(returnValue.sourceHeight, surfaces_Properties[d].sourceHeight);
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();
        layer_dim[c][0] = 85 * (c + 1);
        layer_dim[c][0] = 75 * (c + 1);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0],
                                               layer_dim[c][1]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 12 * (f + 1);
        layer_pos[f][1] = 17 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                ilmSurfaceProperties returnValue;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                ASSERT_EQ(ILM_SUCCESS,
                          ilm_getPropertiesOfSurface(IDs[q], &returnValue));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[s] - layers[0];

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 10 + (b * 35);
        surf_dim[b][1] = 15 + (b * 45);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0] = 85 * (c + 1),
                                               layer_dim[c][1] = 75 * (c + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 12 * (f + 1);
        layer_pos[f][1] = 17 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set source rectangle of layers
    for (int d = 0; d < no_layers; d++)
    {
        layers_Properties[d].sourceX = 79 + (d * 10);
        layers_Properties[d].sourceY = 6438 + (d * 3);
        layers_Properties[d].sourceWidth = 538 + (d * 7);
        layers_Properties[d].sourceHeight = 2 + (d * 2);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerSetSourceRectangle(layers[d],
                                              layers_Properties[d].sourceX,
                                              layers_Properties[d].sourceY,
                                              layers_Properties[d].sourceWidth,
                                              layers_Properties[d].sourceHeight));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Confirm source rectangles of layers
    for (int d = 0; d < no_layers; d++)
    {
        ilmLayerProperties returnValue;
        ASSERT_EQ(ILM_SUCCESS,
                  ilm_getPropertiesOfLayer(layers[d], &returnValue));

        ASSERT_EQ(returnValue.sourceX, layers_Properties[d].sourceX);
        ASSERT_EQ(returnValue.sourceY, layers_Properties[d].sourceY);
        ASSERT_EQ(returnValue.sourceWidth, layers_Properties[d].sourceWidth);
        ASSERT_EQ(returnValue.sourceHeight, layers_Properties[d].sourceHeight);
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[s] - layers[0];
            ilmLayerProperties returnValue;

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfLayer(IDs[s], &returnValue));
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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 8 + (b * 15);
        surf_dim[b][1] = 9 + (b * 35);

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     ILM_PIXELFORMAT_RGBA_8888,
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0] = 25 * (c + 1),
                                               layer_dim[c][1] = 85 * (c + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 14 * (f + 1);
        layer_pos[f][1] = 16 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 12 * (g * 32);
        surf_pos[g][1] = 13 * (g * 48);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Loop round and take layer shots
    for (int h = 0; h < no_layers; h++)
    {
        // make sure the file is not there before
        FILE* f = fopen(outputFile, "r");
        if (f!=NULL){
            fclose(f);
            int result = remove(outputFile);
            ASSERT_EQ(0, result);
        }

        ASSERT_EQ(ILM_SUCCESS, ilm_takeLayerScreenshot(outputFile, layers[h]));

        sleep(1);
        f = fopen(outputFile, "r");
        ASSERT_TRUE(f!=NULL);
        fclose(f);
        remove(outputFile);
    }

    // Loop round and take surface shots
    for (int i = 0; i < no_surfaces; i++)
    {
        const char* outputFile = "/tmp/test.bmp";

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
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[s] - layers[0];
            ilmLayerProperties returnValue;

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfLayer(IDs[s], &returnValue));

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
    for (int b = 0; b < no_surfaces; b++)
    {
        surfaces[b] = getSurface();
        surfaces_dest[b] = surfaces[b];
        surf_dim[b][0] = 18 + (b * 12);
        surf_dim[b][1] = 29 + (b * 37);
        surface_formats[b] = formats[b/8];

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_surfaceCreate((t_ilm_nativehandle)wlSurfaces[b],
                                     surf_dim[b][0],
                                     surf_dim[b][1],
                                     surface_formats[b],
                                     &(surfaces_dest[b])));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

    }

    // Set dimensions of surfaces
    for (int a = 0; a < no_surfaces; a++)
    {
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetDimension(surfaces[a], surf_dim[a]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Create layers
    for (int c = 0; c < no_layers; c++)
    {
        layers[c] = getLayer();

        ASSERT_EQ(ILM_SUCCESS,
                  ilm_layerCreateWithDimension(&(layers[c]),
                                               layer_dim[c][0] = 32 * (c + 1),
                                               layer_dim[c][1] = 72 * (c + 1)));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of layers
    for (int f = 0; f < no_layers; f++)
    {
        layer_pos[f][0] = 14 * (f + 1);
        layer_pos[f][1] = 16 * (f + 1);
        ASSERT_EQ(ILM_SUCCESS, ilm_layerSetPosition(layers[f], layer_pos[f]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Set position of surfaces
    for (int g = 0; g < no_surfaces; g++)
    {
        surf_pos[g][0] = 15 * (g * 38);
        surf_pos[g][1] = 17 * (g * 42);
        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceSetPosition(surfaces_dest[g], surf_pos[g]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());
    }

    // Check pixelmap of surfaces
    for (int h = 0; h < no_surfaces; h++)
    {
        ilmPixelFormat pf;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetPixelformat(surfaces_dest[h], &pf));
        ASSERT_EQ(surface_formats[h], pf);
    }

    // Loop through surfaces and remove
    for (int p = 0; p < no_surfaces; p++)
    {
        t_ilm_int length;
        t_ilm_surface* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_surfaceRemove(surfaces_dest[p]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining surfaces
        ASSERT_EQ(ILM_SUCCESS, ilm_getSurfaceIDs(&length, &IDs));

        // Loop through remaining surfaces and confirm dimensions, positions
        // & Orientations are unchanged
        for (int q = 0; q < length; q++)
        {
            uint index = no_surfaces;
            for (int t = 0; t < no_surfaces; t++)
            {
                if (IDs[q] == surfaces_dest[t]) index = t;
            }

            if ( index != no_surfaces)
            {
                t_ilm_uint dimreturned[2];
                t_ilm_uint posreturned[2];
                ilmPixelFormat pf_returned;
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetDimension(IDs[q], dimreturned));
                EXPECT_EQ(surf_dim[index][0], dimreturned[0]);
                EXPECT_EQ(surf_dim[index][1], dimreturned[1]);
                EXPECT_EQ(ILM_SUCCESS, ilm_surfaceGetPosition(IDs[q], posreturned));
                EXPECT_EQ(surf_pos[index][0], posreturned[0]);
                EXPECT_EQ(surf_pos[index][1], posreturned[1]);
                ASSERT_EQ(ILM_SUCCESS, ilm_surfaceGetPixelformat(IDs[q], &pf_returned));
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
    for (int r = 0; r < no_layers; r++)
    {
        t_ilm_int length;
        t_ilm_layer* IDs;

        ASSERT_EQ(ILM_SUCCESS, ilm_layerRemove(layers[r]));
        ASSERT_EQ(ILM_SUCCESS, ilm_commitChanges());

        // Get remaining layers
        ASSERT_EQ(ILM_SUCCESS, ilm_getLayerIDs(&length, &IDs));

        // Iterate round remaining layers and check positions & dimensions
        for (int s = 0; s < length; s++)
        {
            t_ilm_uint posreturned[2] = {0, 0};
            t_ilm_uint dimreturned[2] = {0, 0};
            t_ilm_int index = IDs[s] - layers[0];
            ilmLayerProperties returnValue;

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetPosition(IDs[s], posreturned));

            EXPECT_EQ(layer_pos[index][0], posreturned[0]);
            EXPECT_EQ(layer_pos[index][1], posreturned[1]);

            EXPECT_EQ(ILM_SUCCESS, ilm_layerGetDimension(IDs[s], dimreturned));

            EXPECT_EQ(layer_dim[index][0], dimreturned[0]);
            EXPECT_EQ(layer_dim[index][1], dimreturned[1]);
            ASSERT_EQ(ILM_SUCCESS,
                      ilm_getPropertiesOfLayer(IDs[s], &returnValue));

        }

        free(IDs);
    }
}





















