#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing this out.

Extract, reproject, resample, and recalculate Landscan data set. 

Created on Thu Mar  5 08:28:07 2020

@author: twillia2
"""

import numpy as np
import subprocess as sp
import rasterio
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path, translate, warp, gdal_options

dpp = Data_Path("/scratch/twillia2/weto/populations/data")
dpc = Data_Path("/projects/rev/data/conus/")

# These come in ESRI grids, we can translate them to geotiffs
sample_src = '~/Downloads/landscan2016_cyprus/cyprus_pop'
sample_dst = dpp.join("test_translation.tif")
translate(sample_src, sample_dst, overwrite=True, format="GTiff")

# We want our geometries to match the exclusion datasets, like this one
template = rasterio.open(dpc.join("_windready_conus.tif"))

# Warp new sample to albers with nearest neighbor first, don't resample
src = dpp.join("rasters/wgs/landscan_2017/landscan_night_2017.tif")
dst = dpp.join("rasters/albers/landscan_near.tif")

# After a few resampling methods, nearest neighbor looks good enough
print(gdal_options("warp"))
extent = list(template.bounds)
srs = template.crs.to_proj4()
warp(src, dst, outputBounds=extent, dstSRS=srs, resampleAlg="near")

# Okay, so we need to resample to 90 meters with mean
src = dpp.join("rasters/albers/landscan_near.tif")
dst = dpp.join("rasters/albers/landscan_near_90m.tif")
warp(src, dst, xRes=90, yRes=90, resampleAlg="average")

# Next, find our scalar and get input on this
original = rasterio.open(src)
transformed = rasterio.open(dst)

# Resolutions
ores = original.res
tres = transformed.res

# Areas
oarea = np.dot(*ores)
tarea = np.dot(*tres)

# Ratios
aratio = tarea / oarea

# Area Ratio overestimates
rratio = tres[0] / ores[0]

# Try this and test the results
src = dst
dst = dpp.join("rasters/albers/landscan_near_90m_count.tif")
sp.call(["gdal_calc.py",
         "-A", src,
         "--calc=(A*" + str(rratio) + ")",
         "--outfile=" + dst,
         "--overwrite"])

# Okay, sum the two up and see how close we are!
src = dpp.join("rasters/wgs/landscan_2017/landscan_night_2017.tif")
dst = dpp.join("rasters/albers/landscan_near_90m_count.tif")

chunks = (1, 5000, 5000)

original = xr.open_rasterio(src, chunks=chunks)[0]
new = xr.open_rasterio(dst, chunks=chunks)[0]

ona = original.attrs["nodatavals"][0]
nna = new.attrs["nodatavals"][0]

with Client():
    
    odata = original.data
    ndata = new.data

    odata[odata == ona] = 0
    ndata[ndata == nna] = 0

    osum = odata.sum().compute()
    nsum = ndata.sum().compute()

print(str(nsum / osum))
