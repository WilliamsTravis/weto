#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add the exculsion layer to our cost code layer to create a coverage map.

Created on Fri Feb  7 12:17:18 2020

@author: twillia2
"""
import numpy as np
import os
import dask.array as da
import xarray as xr
from osgeo import gdal
from weto.functions import Data_Path, warp, to_raster, gmask
from dask.distributed import Client

# data path
dp = Data_Path("/scratch/twillia2/weto/data")

try:
    assert os.path.exists(dp.join("rasters/albers/cost_codes.tif"))
except AssertionError:
    print("The cost code layer has not been generated yet.")

# Check for the exclusion layer
if not os.path.exists(dp.join("rasters/albers/masks/exclusions.tif")):
    src = "/projects/rev/data/conus/_windready_conus.tif"
    dst = "/scratch/twillia2/weto/data/rasters/albers/masks/exclusions.tif"
    template = "/scratch/twillia2/weto/data/rasters/albers/cost_codes.tif"
    warp(src, dst, template)

# Cost Code and exclusion layers
codes = xr.open_rasterio(dp.join("rasters/albers/cost_codes.tif"))[0]
excls = xr.open_rasterio(dp.join("rasters/albers/exclusions.tif"))[0]
dcodes = da.from_array(codes)
dexcls = da.from_array(excls)

# We want a direct mask of dcodes, and an inverse mask of dexcls
mcodes = gmask(dcodes)
mexcls = gmask(dexcls)

# Add these together for the full mask?
coverage = mcodes + mexcls

# And this will have 0s, 1s, and 2s - make it just 1s
coverage = gmask(coverage)

# Now calculate everything
cov = coverage.compute()

# save to raster
reference = gdal.Open(dp.join("rasters/albers/blm_codes_ac.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()
to_raster(cov, dp.join("rasters/albers/masks/coverage.tif"), proj, geom)
