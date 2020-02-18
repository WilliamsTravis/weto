#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a mask of just CONUS cells.

Created on Fri Feb  7 15:42:10 2020

@author: twillia2
"""

import os
import dask.array as da
import xarray as xr
from osgeo import gdal
from weto.functions import Data_Path, to_raster, gmask, warp

# data path
dp = Data_Path("/scratch/twillia2/weto/data")

# The nlcd grid should have every points
src = dp.join("rasters/nlcd_2016_ag.tif")
dst = dp.join("rasters/albers/masks/conus.tif")
template = "/scratch/twillia2/weto/data/rasters/albers/cost_codes.tif"
warp(src, dst, template)
