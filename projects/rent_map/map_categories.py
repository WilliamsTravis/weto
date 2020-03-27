#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a code raster of categories from the original.

Created on Thu Mar 26 12:40:25 2020

@author: twillia2
"""

import numpy as np
import pandas as pd
import rasterio
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path, to_raster

DP = Data_Path("/scratch/twillia2/weto/data")

CATEGORIES = {1: "BLM",
              2: "Tribal", 
              3: "State", 
              4: "Private"}

CHUNKS = {'band': 1, 'x': 5000, 'y': 5000}


def categorize(row):
    if "BLM" in row["type"]:
        return 1
    if "Tribal" in row["type"]:
        return 2
    if "State Land" in row["type"]:
        return 3
    else:
        return 4

    
def main():
    
    # Get lookup table with all of the codes
    lookup = pd.read_csv(DP.join("tables", "conus_cbe_lookup.csv"))
    lookup["category"] = lookup.apply(categorize, axis=1)    

    # Get all codes for each category
    cat_codes = {}
    for cat in CATEGORIES.keys():
        code_vals = lookup["code"][lookup["category"] == cat].values
        cat_codes[cat] = code_vals

    # Open the rasters
    code_path = DP.join("rasters", "albers", "acre", "cost_codes.tif")
    cost_path = DP.join("rasters", "albers", "acre", "rent_map.tif")
    template = rasterio.open(code_path)
    codes = xr.open_rasterio(code_path, chunks=CHUNKS)[0].data
    costs = xr.open_rasterio(cost_path, chunks=CHUNKS)[0].data 
    new_costs = costs.copy()
    new_codes = codes.copy()

    # Okay, change all of the values
    for cat, code_vals in cat_codes.items():
        minc = np.nanmin(code_vals)
        maxc = np.nanmax(code_vals)
#        new_costs[(codes >= minc) & (codes <= maxc) & (costs > 0)] = cat
        new_codes[(codes >= minc) & (codes <= maxc)] = cat
        
    # Now compute the results
    with Client():
#        new_costs = new_costs.compute()
        new_codes = new_codes.compute()

    # Save
    proj = template.crs.to_wkt()
    geom = template.transform.to_gdal()
    new_cost_path = DP.join("rasters", "albers", "acre", "cost_cats.tif")
    new_code_path = DP.join("rasters", "albers", "acre", "code_cats.tif")
    to_raster(new_costs, new_cost_path, proj, geom, dtype="int16",
              compress="LZW")
    to_raster(new_codes, new_code_path, proj, geom, dtype="int16",
              compress="LZW")

if __name__ == "__main__":
    main()
