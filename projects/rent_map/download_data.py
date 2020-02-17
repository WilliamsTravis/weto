#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Downloading and transforming data.

Created on Mon Feb  3 12:35:54 2020

@author: travis
"""
import xarray as xr
from weto.gdalmethods import to_raster

# NLCD
nlcd_path = ("https://s3-us-west-2.amazonaws.com/mrlc/NLCD_2016_Land_Cover"
             "_L48_20190424.zip")

# BLM
blm_path = ""

# Tribal Land
tribal_path = ""

# State Land
state_path =  ""

# Counte Boundaries
cnty_path = ""
