#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate coverage by census blocks.

Created on Tue Mar 24 10:23:16 2020

@author: twillia2
"""

import geopandas as gpd
from coverage import codes

# Get a shapefile of the census blocks
regions = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2018/shp/"
                        "cb_2018_us_region_5m.zip")
divisions = gpd.read_file("https://www2.census.gov/geo/tiger/GENZ2018/shp/"
                          "cb_2018_us_division_5m.zip")

# Get the codes for each category
codes = codes()

