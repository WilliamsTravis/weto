#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
We need a mask for just full exclusions. 1s for full exclusions 0s for
everything else.

Created on Fri Feb 14 09:18:28 2020

@author: twillia2
"""
import dask.array as da
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path, to_raster, warp

# Data Paths
# dp = Data_Path("~/Box/WETO 1.2/data")
dp = Data_Path("/scratch/twillia2/weto/data")
exl_path = dp.join("rasters/core_exclusions_raster/alopez_core_exclusions.tif")

# Best chunk size?
chunks = {'band': 1, 'x': 5000, 'y': 5000}

# Using the 'core exclusions raster' that antyhony put together. Check this.
xexcl = xr.open_rasterio(exl_path, chunks=chunks)[0]
excl = xexcl.data

# We need to reverse this
excl[da.isnan(excl)] = 0
excl = (excl - 1) * -1

# And let's make exclusion values 9999
excl[excl == 1] = 9999

# Compute
with Client():
    excl = excl.compute()

# save to raster
to_raster(excl, dp.join("rasters/inverse_exclusions.tif"), template=exl_path,
          compress="deflate")

# warp to acre grid in north american albers equal area conic
res = 63.614907234075254
warp(dp.join("rasters/inverse_exclusions.tif"),
     dp.join("rasters/albers/acre/inverse_exclusions.tif"),
     xRes=res,
     yRes=res,
     overwrite=True)