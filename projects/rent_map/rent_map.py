#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map dollar values to the codes in cost_codes.tif

Created on Fri Feb 21 10:02:49 2020

@author: twillia2
"""
import numpy as np
import os
import pandas as pd
import subprocess as sp
from gdalmethods import Data_Path, Map_Values
from glob import glob

# Set Data Path
dp = Data_Path("/scratch/twillia2/weto/data")

# Translator for the lookup table
def fixit(x):
    x = x.replace("*", "").replace("\n", "").replace("$", "").replace(",", "")
    if x == " NULL " or  x == "  NULL  ":
        x = 0
    return x

# Create the lookup table
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
lookup.columns = ['code', 'type', 'dollar_ac']
lookup["dollar_ac"] = lookup["dollar_ac"].apply(fixit)
lookup["dollar_ac"] = lookup["dollar_ac"].astype(float)

# Tile
code_path = dp.join("rasters/albers/acre/cost_codes.tif")
code_tile_folder = dp.join("rasters/albers/acre/code_tiles")
os.makedirs(code_tile_folder, exist_ok=True)
sp.call(["gdal_retile.py",
         "-ps", "5000", "5000",
         "-targetDir", code_tile_folder,
         code_path])
code_tiles = glob(os.path.join(code_tile_folder, "*tif"))

# Map
cost_tile_folder = dp.join("rasters/albers/acre/cost_tiles")
val_dict = dict(zip(lookup["code"], lookup["dollar_ac"]))
mv = Map_Values(val_dict)
cost_tiles = mv.map_files(code_tiles, cost_tile_folder, ncpu=15)

# Merge
cost_file = dp.join("rasters/albers/acre/rent_map.tif")
call =  ["gdal_merge.py", "-o", cost_file] + cost_tiles
sp.call(call)

# Done.
