#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Taked the masked code grids, stack them, and create a single code grid.

There is obviously a better way.

Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""
import dask.array as da
import xarray as xr
from osgeo import gdal
from weto.functions import Data_Path, to_raster

# data path
dp = Data_Path("/scratch/twillia2/weto/data")

# Best chunk size?
chunks = {'band': 1, 'x': 5000, 'y': 5000}

# Get each array, remove band dimension
blm = xr.open_rasterio(dp.join("rasters/albers/blm_codes_ac.tif"))[0]  # <----- This needs to be masked by exclusions.
tribes = xr.open_rasterio(dp.join("rasters/albers/masks/masked_tribes.tif"))[0]
state = xr.open_rasterio(dp.join("rasters/albers/masks/masked_state.tif"))[0]
nlcd = xr.open_rasterio(dp.join("rasters/albers/masks/masked_nlcd.tif"))[0]
excl = xr.open_rasterio(dp.join("data/rasters/albers/acre/exclusions.tif"))[0] # Needs to have 1s for full exclusions and 0s for all else.

# Convert to dask arrays ?
dnlcd = da.from_array(nlcd)
dblm = da.from_array(blm)
dtribes = da.from_array(tribes)
dstate = da.from_array(state)

# ...Create an ndarray out of the originals, and take the max again
composite = da.stack([excl, dblm, dtribes, dstate, dnlcd], axis=0).max(axis=0)
final = composite.compute()  # memory explosion?

# Save to file
reference = gdal.Open(dp.join("rasters/albers/blm_codes_ac.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()
to_raster(final, dp.join("rasters/albers/cost_codes.tif"), proj, geom)

# Done?
