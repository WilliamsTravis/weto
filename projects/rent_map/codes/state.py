#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rasterize State Land Codes.

Created on Mon Feb  3 10:35:24 2020

@author: twillia2
"""

import geopandas as gpd
import os
import pandas as pd
from osgeo import gdal
from gdalmethods import Data_Path, rasterize, reproject_polygon

# Data Path
dp = Data_Path("/scratch/twillia2/weto/data")

# Lookup table
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
lookup.columns = ['code', 'type', 'dollar_ac']

# Get just the state names
state_lu = lookup[["type", "code"]][lookup["type"].str.contains("State Land")]
get_state = lambda x: x[:x.index("State Land") - 1]
state_lu["state_nm"] = state_lu["type"].apply(get_state)
state_lu["state_nm"] = state_lu["state_nm"].astype(str)

# First add the codes from the lookup table to the shapefiles and rewrite
state = gpd.read_file(dp.join("shapefiles/USA/conus_padus_state.shp"))

# Join state_lu with state
state = state.merge(state_lu, on="state_nm", how="left")
state = state[["gid", "code", "geometry"]]

# Reproject
state_shp = dp.join("shapefiles/USA/conus_padus_codes.shp")
state_shp_albers = dp.join("shapefiles/USA/albers/conus_padus_codes.shp")
t_srs = ("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 "
         "+y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
state.to_file(state_shp)
reproject_polygon(src=state_shp, dst=state_shp_albers, t_srs=t_srs)

# get the target geometry
template = dp.join("rasters/albers/acre/nlcd_ag.tif")

# Rasterize the code field to the finest resolution we'll use
state_tif = dp.join("rasters/albers/acre/state_codes.tif")
if not os.path.exists(state_tif):
    rasterize(src=state_shp_albers,
              dst=state_tif,
              attribute="code",
              template_path=template,
              dtype=gdal.GDT_Float32,
              overwrite=True)
