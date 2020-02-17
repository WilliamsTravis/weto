#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create masks for the ranking procedure. Separated here for memory purposes

Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""
import dask.array as da
import numpy as np
import os
from osgeo import gdal
from weto.functions import Data_Path, to_raster
import xarray as xr


# data path
dp = Data_Path("/scratch/twillia2/weto/data")

# Get each array, remove band dimension
blm = xr.open_rasterio(dp.join("rasters/albers/blm_codes_ac.tif"))[0]
tribes = xr.open_rasterio(dp.join("rasters/albers/tribal_codes_ac.tif"))[0]
state = xr.open_rasterio(dp.join("rasters/albers/state_land_codes_ac.tif"))[0]

# Make a mask of each higher priority layer
bmask = gmask(blm)
tmask = gmask(tribes)
smask = gmask(state)

# Make two composite masks, the max will bring 1s to the front
mask1 = da.stack([bmask, tmask, smask], axis=0).max(axis=0)
mask2 = da.stack([bmask, tmask], axis=0).max(axis=0)
mask3 = bmask

# Now we actually want the inverse of these (1 -> 0 -> 0 ; 0 -> -1 -> 1)
mask1 = (mask1 - 1) * -1
mask2 = (mask2 - 1) * -1
mask3 = (mask3 - 1) * -1

# Compute each of these and save to file
os.makedirs(dp.join("rasters/albers/masks/"), exist_ok=True)
reference = gdal.Open(dp.join("rasters/albers/blm_codes_ac.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()

mask1 = mask1.compute().value.data
to_raster(mask1, dp.join("rasters/albers/masks/mask1.tif"), proj, geom)
del mask1

mask2 = mask2.compute().value.data
to_raster(mask2, dp.join("rasters/albers/masks/mask2.tif"), proj, geom)
del mask2

mask3 = mask3.compute()
to_raster(mask3, dp.join("rasters/albers/masks/mask3.tif"), proj, geom)
del mask3
