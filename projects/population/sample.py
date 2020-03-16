#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Associate landscan population estimates with reV supply curve points.

Created on Thu Feb 27 08:56:48 2020

@author: twillia2
"""
import geopandas as gpd
import numpy as np
import pandas as pd
import subprocess as sp
import matplotlib.pyplot as plt
import os
import rasterio

from gdalmethods import Data_Path, to_geo, warp

# Data Paths
DP = Data_Path("/scratch/twillia2/weto/populations")
DPA = Data_Path(DP.join("albers"))
DPSC = Data_Path("/projects/rev/new_projects/ipm_wind/outputs")
DPTMP = Data_Path("/projects/rev/data/conus")
PROJ = ("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 "
        "+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Reproject Supply Curve Table
if not os.path.exists(DPA.join("sc_gids_singles.gpkg")):
    sc = pd.read_csv(DPSC.join("outputs_sc.csv"))
    sc = sc[["sc_gid", "latitude", "longitude"]]
    sc = sc.drop_duplicates(keep="first", subset=["latitude", "longitude"])
    sc = sc.reset_index(drop=True)
    sc = to_geo(sc, loncol="longitude", latcol="latitude", epsg=4326)
    sc.to_file(DP.join("sc_gids_singles.gpkg"), driver="GPKG")
    src = DP.join("sc_gids_singles.gpkg")
    dst = DPA.join("sc_gids_singles.gpkg")
    sp.call(["ogr2ogr", dst, src, "-t_srs", PROJ])

# Reproject population, keep an eye on the values
if not os.path.exists(DPA.join("landscan_night_2017.tif")):
    src = DP.join("landscan_2017",  "landscan_night_2017.tif")
    dst = DPA.join("landscan_night_2017.tif")
    warp(src, dst, dstSRS=PROJ)  # <------------------------------------------- Default nearest neighbors

# All of our point data sets have irregular spacings! This might need to be projected first
sc = gpd.read_file(DPA.join("sc_gids_singles.gpkg"))
sc["x"] = sc["geometry"].apply(lambda x: x.x)
sc["y"] = sc["geometry"].apply(lambda x: x.y)
ydiffs = np.diff(np.unique(sc["y"]))
plt.hist(ydiffs, bins=1000)
distance = 5700

# Circle buffer
sc_buffer = sc.copy()
sc_buffer["geometry"] = sc_buffer.buffer(distance / 2, cap_style=1)
sc_buffer = sc_buffer[["sc_gid", "geometry"]]
sc_buffer.crs = PROJ
sc_buffer.to_file(DPA.join("sc_circle_single_buffer.gpkg"), driver="GPKG")

# Square Buffer
sc_buffer["geometry"] = sc_buffer.envelope
sc_buffer.to_file(DPA.join("sc_square_single_buffer.gpkg"), driver="GPKG")

# Sum up zonal stats
from rasterstats import zonal_stats
poppath = DPA.join("landscan_night_2017.tif")
p = rasterio.open(poppath)
zonepath = DPA.join("sc_circle_buffer.gpkg")
out = zonal_stats(zonepath, poppath, stats="sum", nodata=p.nodata)
