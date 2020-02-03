#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GDAL Command holder.

Created on Mon Feb  3 09:24:36 2020

@author: twillia2
"""

# To warp .000326 degree wgs raster to albers
gdalwarp -s_srs EPSG:4326 -t_srs "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" /scratch/twillia2/weto/data/rasters/nlcd_index.tif /scratch/twillia2/weto/data/rasters/albers/nlcd_index.tif

# To resample to an acre
63.614907234075254
gdalwarp -tr 63.614907234075254 63.614907234075254 in.tif out.tif