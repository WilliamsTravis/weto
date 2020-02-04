#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GDAL Command holder.

Created on Mon Feb  3 09:24:36 2020

@author: twillia2
"""

# To warp .000326 degree wgs raster to albers at 30m
gdalwarp -tr 30 30 -s_srs EPSG:4326 -t_srs "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" in.tif out.tif

# To resample to an acre
gdalwarp -tr 63.614907234075254 63.614907234075254 in.tif out.tif