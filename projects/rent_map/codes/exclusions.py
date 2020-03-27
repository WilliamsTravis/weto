# -*- coding: utf-8 -*-
"""
We need a mask for just full exclusions. 1s for full exclusions 0s for
everything else.

Created on Fri Feb 14 09:18:28 2020

@author: twillia2
"""

import dask.array as da
import numpy as np
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path, to_raster, warp

# Data Paths
# dp = Data_Path("~/Box/WETO 1.2/data")
DP = Data_Path("/scratch/twillia2/weto/data")
EXL_PATH = DP.join("rasters", "core_exclusions_raster",
                   "alopez_core_exclusions.tif")
ROAD_PATH = DP.join("rasters", "core_exclusions_raster", "conus_roads.tif")
CONUS_PATH = DP.join("rasters", "albers", "90m", "conus.tif")
RAIL_PATH =  DP.join("rasters", "core_exclusions_raster", "conus_rail.tif")

# Best chunk size?
CHUNKS = {'band': 1, 'x': 5000, 'y': 5000}


def build_exclusions():
    """ Build exclusion file for the WETO rent map."""

    # Using the 'core exclusions raster' that antyhony put together.
    excl = xr.open_rasterio(EXL_PATH, chunks=CHUNKS)
    roads = xr.open_rasterio(ROAD_PATH, chunks=CHUNKS)
    conus = xr.open_rasterio(CONUS_PATH, chunks=CHUNKS)
    rails = xr.open_rasterio(RAIL_PATH, chunks=CHUNKS)

    # Different na values every time :/
    roadsna = roads.attrs["nodatavals"][0]
    conusna = conus.attrs["nodatavals"][0]
    railsna = rails.attrs["nodatavals"][0]

    # Get just the data arrays
    excl = excl[0].data
    roads = roads[0].data
    conus = conus[0].data
    rails = rails[0].data

    # Set nodata values to 0
    excl[da.isnan(excl)] = 0
    roads[roads == roadsna] = 0
    conus[conus == conusna] = 0
    rails[rails == railsna] = 0

    # We need to reverse the original exclusions
    excl = (excl - 1) * -1

    # Combine roads
    excl = da.stack([excl, roads, rails], axis=0).max(axis=0)

    # And let's make exclusion values 9999 since 1 will be a code
    excl[excl == 1] = 9999

    # And cut out just CONUS for mapping
    excl = excl * conus

    # Compute
    print("Combining exclusion layers...")
    with Client():
        excl = excl.compute()

    # save to raster
    print("Saving to 90 meter reV grid...")
    to_raster(excl, DP.join("rasters", "rent_exclusions.tif"),
              template=EXL_PATH, compress="deflate")

    # warp to acre grid in north american albers equal area conic
    print("Warping to acre grid...")
    res = 63.614907234075254
    warp(DP.join("rasters", "rent_exclusions.tif"),
         DP.join("rasters", "albers", "acre", "rent_exclusions.tif"),
         xRes=res,
         yRes=res,
         overwrite=True)

    print("Done.")

if __name__ == "__main__":
    build_exclusions()
