#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create masks for the ranking procedure. Separated here for memory purposes.

Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""
import dask.array as da
import numpy as np
import os
import xarray as xr
from dask.distributed import Client
from gdalmethods import Data_Path, to_raster
from osgeo import gdal

# Move to functions
@da.as_gufunc(signature="(i)->(i)", output_dtypes=int, vectorize=True,
              allow_rechunk=True)
def gmask(array, navalues=[-9999., 0.0]):
    """Create a mask of 1's and 0 out of an array using non values"""
    # So multiple values can be masked out
    if not isinstance(navalues, (list, tuple, np.ndarray)):
        navalues = [navalues]

    # Create mask
    array = array.copy()
    array[np.isin(array, navalues)] = np.nan
    array = array * 0 + 1
    array[np.isnan(array)] = 0

    return array

# data paths
dp = Data_Path("/scratch/twillia2/weto/data")
mask_dp = Data_Path(dp.join("rasters/albers/acre/masks"))

# Set the chunk sizes
chunks = {'band': 1, 'x': 5000, 'y': 5000}

# Reference geometry
reference = gdal.Open(dp.join("rasters/albers/acre/blm_codes.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()

# Code paths
nlcd_path = dp.join("rasters/albers/acre/nlcd_codes.tif")
blm_path = dp.join("rasters/albers/acre/blm_codes.tif")
tribal_path = dp.join("rasters/albers/acre/tribal_codes.tif")
state_path = dp.join("rasters/albers/acre/state_codes.tif")
excl_path = dp.join("rasters/albers/acre/inverse_exclusions.tif")
tobe_conus_path = dp.join("rasters/albers/acre/nlcd.tif")

# Get each array, remove band dimension
nlcd = xr.open_rasterio(nlcd_path, chunks=chunks)[0].data
blm = xr.open_rasterio(blm_path, chunks=chunks)[0].data
tribes = xr.open_rasterio(tribal_path, chunks=chunks)[0].data
state = xr.open_rasterio(state_path, chunks=chunks)[0].data
excl = xr.open_rasterio(excl_path, chunks=chunks)[0].data

# I need to fix this...nlcd's shape is off by one cell from my tiling/merging
if nlcd.shape != blm.shape:
    # One too many on top
    nlcd = nlcd[1:, :]

    # Now we need two on bottom
    x = da.from_array(np.zeros((1, nlcd.shape[1])))
    nlcd = np.vstack((nlcd, x, x))

    # and one extra on the side
    y = da.from_array(np.zeros((nlcd.shape[0], 1)))
    nlcd = np.hstack((nlcd, y))
    with Client():
        nlcd = nlcd.compute()
    to_raster(nlcd, nlcd_path, proj, geom)

# Make a mask of each higher priority layer
emask = gmask(excl)
bmask = gmask(blm)
tmask = gmask(tribes)
smask = gmask(state)

# Make a conus mask with 1s for in and nans for out
if not os.path.exists(mask_dp.join("conus.tif")):
    conus = xr.open_rasterio(tobe_conus_path, chunks=chunks)[0].data
    conus[conus > 0] = 1
    conus[conus <= 0] = np.nan
    with Client():
        conus = conus.compute()
        print("Saving Conus")
        to_raster(conus, mask_dp.join("conus.tif"), proj, geom)
        del conus

# Make four composite masks, the max will bring 1s to the front
mask1 = da.stack([emask, bmask, tmask, smask], axis=0).max(axis=0)
mask2 = da.stack([emask, bmask, tmask], axis=0).max(axis=0)
mask3 = da.stack([emask, bmask], axis=0).max(axis=0)
mask4 = emask

# Now we actually want the inverse of these (1 -> 0 -> 0 ; 0 -> -1 -> 1)
mask1 = (mask1 - 1) * -1
mask2 = (mask2 - 1) * -1
mask3 = (mask3 - 1) * -1
mask4 = (mask4 - 1) * -1

# We'll need to setup interim files, I can't handle the memory useage
layers = {"nlcd": nlcd,
          "state": state,
          "tribes": tribes,
          "blm": blm}
masked_paths = {"nlcd": mask_dp.join("nlcd.tif"),
                "state": mask_dp.join("state.tif"),
                "tribes": mask_dp.join("tribes.tif"),
                "blm": mask_dp.join("blm.tif")}
masks =  {"nlcd": mask1,
          "state": mask2,
          "tribes": mask3,
          "blm": mask4}

# Loop through, the lowest priority layers are masked by the largest mask <---- This is obviously not the best way to do this
for key, layer in layers.items():
    layer_path = masked_paths[key]
    mask = masks[key]
    mlayer = layer * mask
    with Client():
        masked_layer = mlayer.compute()
    to_raster(masked_layer, layer_path, proj, geom)
    del masked_layer

# Too much at once, read the precomputd masked layers from file
mblm = xr.open_rasterio(masked_paths["blm"], chunks=chunks)[0].data
mtribes = xr.open_rasterio(masked_paths["tribes"], chunks=chunks)[0].data
mstate = xr.open_rasterio(masked_paths["state"], chunks=chunks)[0].data
mnlcd = xr.open_rasterio(masked_paths["nlcd"], chunks=chunks)[0].data
mconus = xr.open_rasterio(mask_dp.join("conus.tif"), chunks=chunks)[0].data

with Client():
    layers = [excl, mblm, mtribes, mstate, mnlcd]
    composite = da.stack(layers, axis=0).max(axis=0)
    composite = composite * mconus
    final = composite.compute()  # memory explosion?

# Save to file
reference = gdal.Open(dp.join("rasters/albers/acre/blm_codes.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()
to_raster(final, dp.join("rasters/albers/acre/cost_codes.tif"), proj, geom)
