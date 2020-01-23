#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create data layers for calculating costs in the LandBOSSE Wind Production
model.

Notes:
    - Match geometry of WTK resource layers. So, if LandBOSSE uses the HDF5
      files like reV does, this is as simple as associating data with the
      WTK grid IDs...?
    - The task calls for a stacked cost surface raster to start, then it
      suggests an HDF5 to incorporate into reV. Is the stacked raster
      necessary?
    - The data sets to transform will be stored somewhere in:
        ~/Box/WETO 1.2/data


Created on Wed Jan  8 11:14:54 2020

@author: twillia2
"""
import os
from geofeather import from_geofeather, to_geofeather
import geopandas as gpd
import pandas as pd
from glob import glob
from functions import shp_to_h5

# I might share this on github
os.chdir(os.path.expanduser("~/github/weto"))

# Expand pandas printouts
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# A sample shapefile layer
shp = "/Users/twillia2/Box/WETO 1.2/data/BLM/fedlands_blm.shp"
trgt = "data/tables/wtk_fedland.h5"

# Shapefile to HDF5 data set
gdf = gpd.read_file(shp)
shp_to_h5(gdf, trgt, attribute="raster_val", mode="a")




