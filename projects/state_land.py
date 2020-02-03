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
import rasterio
from osgeo import gdal
from weto.functions import Data_Path, rasterize


# Data Path
dp = Data_Path("/scratch/twillia2/weto/data")

# Lookup table
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
lookup.columns = ['code', 'type',' dollar_ac']

# First add the codes from the lookup table to the shapefiles and rewrite
state = gpd.read_file(dp.join("shapefiles/USA/tl_2016_us_aiannh.shp"))

# This one only has one value
codes = lookup["code"][lookup["type"] == "Tribal Land"].values[0]

# Add this code into the shapefile and rewrite
tribal["code"] = code
tribal.to_file(dp.join("shapefiles/tribal/tl_2016_us_aiannh.shp"))

# get the target geometry
nlcd = rasterio.open(dp.join("rasters/nlcd_2016_ag.tif"))
transform = nlcd.get_transform()
height = nlcd.height
width = nlcd.width

# Rasterize the code field to the finest resolution we'll use
if not os.path.exists(dp.join("rasters/tribal_codes.tif")):
    rasterize(src=dp.join("shapefiles/tribal/tl_2016_us_aiannh.shp"),
              dst=dp.join("rasters/tribal_codes.tif"),
              attribute="code",
              transform=transform,
              height=height,
              width=width,
              epsg=4326,
              dtype=gdal.GDT_Float32,
              overwrdite=True)
