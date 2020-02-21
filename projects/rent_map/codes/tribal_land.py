#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rasterize tribal land.

Created on Mon Feb  3 09:30:06 2020

@author: twillia2
"""

import geopandas as gpd
import pandas as pd
from osgeo import gdal
from gdalmethods import Data_Path, rasterize, reproject_polygon


# Data Path
dp = Data_Path("/scratch/twillia2/weto/data")

# Lookup table
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
lookup.columns = ['code', 'type', 'dollar_ac']

# First add the codes from the lookup table to the shapefiles and rewrite
tribal = gpd.read_file(dp.join("shapefiles/tribal/tl_2016_us_aiannh.shp"))

# This one only has one value
code = lookup["code"][lookup["type"] == "Tribal Land"].values[0]

# Add this code into the shapefile and reproject
tribal["code"] = code
shp_path = dp.join("shapefiles/tribal/tribal_codes.shp")
shp_path_albers = dp.join("shapefiles/tribal/albers/tribal_codes.shp")
t_srs = ("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 "
         "+y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") 
tribal.to_file(shp_path)
reproject_polygon(shp_path, shp_path_albers, t_srs=t_srs)

# get the target geometry
template = dp.join("rasters/albers/acre/nlcd_ag.tif")

# Rasterize the code field to the finest resolution we'll use
rasterize(src=shp_path_albers,
          dst=dp.join("rasters/albers/acre/tribal_codes.tif"),
          attribute="code",
          template_path=template,
          dtype=gdal.GDT_Float32,
          overwrite=True)
