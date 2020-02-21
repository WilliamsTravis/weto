#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply the masks from rank_mask.py to the code grids.

Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""
import dask.array as da
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
nlcd = xr.open_rasterio(dp.join("rasters/albers/nlcd_codes_ac.tif"))[0]

# convert to dask arrays
dnlcd = da.from_array(nlcd)
dblm = da.from_array(blm)
dtribes = da.from_array(tribes)
dstate = da.from_array(state)

# Read masks
mask1 = xr.open_rasterio(dp.join("rasters/albers/masks/mask1.tif"))[0]
mask2 = xr.open_rasterio(dp.join("rasters/albers/masks/mask2.tif"))[0]
mask3 = xr.open_rasterio(dp.join("rasters/albers/masks/mask3.tif"))[0]

# convert to dask arrays
dmask1 = da.from_array(mask1)
dmask2 = da.from_array(mask2)
dmask3 = da.from_array(mask3)

# Then we can remove all overlapping values in the order required...?
mnlcd = dnlcd * dmask1
mstate = dstate * dmask2
mtribes = dtribes * dmask3

# Compute each of these and save to file
os.makedirs(dp.join("rasters/albers/masks/"), exist_ok=True)
reference = gdal.Open(dp.join("rasters/albers/blm_codes_ac.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()

mnlcd = mnlcd.compute()
to_raster(mnlcd, dp.join("rasters/albers/masks/masked_nlcd.tif"), proj, geom)
del mnlcd

mstate = mstate.compute()
to_raster(mstate, dp.join("rasters/albers/masks/masked_state.tif"), proj, geom)
del mstate

mtribes = mtribes.compute()
to_raster(mtribes, dp.join("rasters/albers/masks/masked_tribes.tif"), proj, geom)
del mtribes
