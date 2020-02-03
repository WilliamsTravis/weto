#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rasterize BLM zone and state land lookup codes.

1) Get the BLM zone shapefile/
2) Join the code value to each zone from our lookup.
3) Use the geometry of the WGS NLCD codes GeoTiff to rasterize
4) Reproject into Albers Equal Area Conic (EPSG 102008)
5) I think if we get the 30m raster, the hard part will be done and we can
   resample to whatever resolution later.

Created on Mon Feb  3 09:30:06 2020

@author: twillia2
"""

import os
import pandas as pd
import rasterio
from osgeo import gdal
from weto.functions import Data_Path, rasterize

import geopandas as gpd

# Data Path
dp = Data_Path("/scratch/twillia2/weto/data")

# Lookup table
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
lookup.columns = ['code', 'type',' dollar_ac']

# First add the codes from the lookup table to the shapefiles and rewrite
blm = gpd.read_file(dp.join("shapefiles/BLM/conus_fedland_blm_county.shp"))

# Which field corresponds to zones?

# Join lookup table to the the blm gdf using the Zone field


# get the target geometry
nlcd = rasterio.open(dp.join("rasters/nlcd_2016_ag.tif"))
transform = nlcd.get_transform()
height = nlcd.height
width = nlcd.width

# Rasterize the code field to the finest resolution we'll use
if not os.path.exists(dp.join("rasters/county_gids.tif")):
    rasterize(src=dp.join("shapefiles/USA/conus_fedland_blm_county.shp"),
              dst=dp.join("rasters/blm_codes.tif"),
              attribute="code",
              transform=transform,
              height=height,
              width=width,
              epsg=4326,
              dtype=gdal.GDT_Float32,
              overwrite=True)
