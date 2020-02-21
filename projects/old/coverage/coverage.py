#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a coverage layer with exclusions.

1) full and partial exclusions...don't count partial exclusions?
2) NA's for non-CONUS pixels
3) 1 for codes and for full exclusions4) 0 for no codes and for partial exclusions

Created on Sun Feb  9 22:33:54 2020

@author: twillia2
"""
import dask.array as da
import numpy as np
import os
import pandas as pd
import subprocess as sp
import xarray as xr
from gdalmethods import Data_Path, warp, to_raster
from dask.distributed import Client
from osgeo import gdal

# Data path
dp = Data_Path("~/github/weto/data")

# paths
conus_path = dp.join("rasters/albers/acre/conus.tif")
codes_path = dp.join("rasters/albers/acre/cost_codes_ac.tif")
excl_path = dp.join("rasters/albers/acre/exclusions_ac.tif")
incl_path = dp.join("rasters/albers/acre/inclusion_ac.tif")
coverage_path = dp.join("rasters/albers/acre/coverage.tif")

# Best chunk size?
chunks = {'band': 1, 'x': 5000, 'y': 5000}

# Create a CONUS mask
if not os.path.exists(conus_path):

    # Use the nlcd_2016 raster, that covers everything.
    src = dp.join("rasters/nlcd_2016.tif")
    dst = dp.join("rasters/albers/acre/nlcd_2016.tif")

    # Create the new nlcd file
    warp(src, dst, template=codes_path, dtype="Int32", src_nodata=0,
         dst_nodata=-9999, overwrite=True)

    # Create field of 1s and 0s using gdal_calc.py
    sp.call(["gdal_calc.py",
             "-A", dst,
             "--outfile", conus_path,
             "--type=", gdal.GDT_Float32,
             "--calc=((A*0)+1)",
             "--NoDataValue=-9999."])

# Create an inclusion raster (opposite of the exclusion with no partials)
if not os.path.exists(incl_path):

    # Read in exclusions layer
    excl = xr.open_rasterio(excl_path)[0].data

    # Remove 0.5 values
    excl[excl == 0.5] = 1
    # excl[excl == 0.5] = 0  # To consider partial exclusions as well

    # Generate the inverse of the resulting exlusions layer
    incl = (excl - 1) * -1

    # Mask by conus
    conus = xr.open_rasterio(conus_path)[0]
    conus = da.from_array(conus)
    incl = incl * conus

    # Save this to file
    incl = incl.compute()
    to_raster(incl, incl_path, template=codes_path, compress="deflate")

if not os.path.exists(coverage_path):
    # Read in our three rasters - xarray seems to expect at least one time dim
    codes = xr.open_rasterio(codes_path, chunks=chunks)[0].data
    incl = xr.open_rasterio(incl_path, chunks=chunks)[0].data
    conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data

    # Now create the overall coverage layer with incl and codes
    conus[conus == -9999.] = np.nan
    coverage = da.stack([codes,  incl], axis=0).max(axis=0)
    coverage = coverage * conus
    with Client():
        coverage = coverage.compute()

    to_raster(coverage, coverage_path, template=codes_path, tiled=False,
              compress="deflate")


if not os.path.exists(dp.join("tables/conus_cbe_lookup.csv")):
    # Fix it...
    print("no lookup table")

# Now we'll need each of these

coverage = xr.open_rasterio(coverage_path, chunks=chunks)[0].data
conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data
incl = xr.open_rasterio(incl_path, chunks=chunks)[0].data
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))

# Codes
blm = lookup[lookup["type"].str.contains("BLM Zone")]
state = lookup[lookup["type"].str.contains("State Land")]
private = lookup[(~lookup["type"].str.contains("BLM Zone")) &
                 (~lookup["type"].str.contains("Tribal Land")) &
                 (~lookup["type"].str.contains("State Land"))]
blm_codes = blm["code"].values
state_codes = state["code"].values
private_codes = private["code"].values

# What is the plan here...how much of each category is represented, or how
# how much of what is represented is in each category? Also, we need to be
# careful not to include exclusions when counting each of the categories.
# In fact, if I was to go back and rank exclusions as the top priority
# in the first steps this would be taken care of and we could skip most of the
# mess above!

# Call these one at a time to avoid overly complicated dask graphs
with Client():
    # Total number of cells within conus
    conus[conus == -9999.] = 0
    ntotal = da.count_nonzero(conus).compute()

# Number of sites covered by codes or full exclusions
with Client():
    coverage[da.isnan(coverage)] = 0
    ncovered = da.count_nonzero(coverage).compute()

# BLM Coverage - we'll need the original rasters of each of these layers
with Client():
    # blm_total = ...`
    blm_covered = coverage[da.isin(coverage, blm_codes)]
    nblm = da.count_nonzero(blm_covered).compute()

# State Percentage
with Client():
    # state_total = ...
    state_covered = coverage[da.isin(coverage, state_codes)]
    nstate = da.count_nonzero(state_covered).compute()

# Private Percentage
with Client():
    # private_total = ...
    private_covered = coverage[da.isin(coverage, private_codes)]
    nprivate = da.count_nonzero(private_covered).compute()

with Client():
    # tribal total = ...
    tribal_covered = coverage[coverage == 16]
    ntribal = da.count_nonzero(tribal_covered).compute()

with Client():
    # exclusions
    # incl_covered = coverage[coverage == 1]
    incl_covered = incl[incl == 1]
    nincl = da.count_nonzero(incl_covered).compute()

# Total coverage
pcoverage = ncovered / ntotal
pblm = nblm / ntotal
pstate = nstate / ntotal
pprivate = nprivate / ntotal
ptribal = ntribal / ntotal
pexcl = nincl / ntotal
