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
from functions import to_geo, rasterize, reproject_polygon, reproject_point

# I might share this on github
os.chdir(os.path.expanduser("~/github/weto"))

# Expand pandas printouts
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


# 1) Create Grids of WTK GIDS
# Point file saved as geofeather for faster I/O
if not os.path.exists("data/shapefiles/wtk_points.feather"):
    wtkpnts = pd.read_csv("data/tables/wind_conus_v1_coords.csv")
    wtkpnts["gid"] = wtkpnts.index
    wtk = to_geo(wtkpnts, loncol="lon", latcol="lat", epsg=4326)
    to_geofeather(wtk, "data/shapefiles/wtk_points.feather")
else:
    wtk = from_geofeather("data/shapefiles/wtk_points.feather")

# Write to shapefile or geopackage
if not os.path.exists("data/shapefiles/wtk_points.shp"):
    wtk.to_file("data/shapefiles/wtk_points.shp")

# Reproject to NA Alber Equal Area Conic <------------------------------------- What do we usually use?
tproj = ("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 " +
         "+y_0=0 +datum=NAD83 +units=m +no_defs")
if not os.path.exists("data/shapefiles/wtk_points_albers.shp"):
    reproject_point(src="data/shapefiles/wtk_points.shp",
                    dst="data/shapefiles/wtk_points_albers.shp",
                    tproj=tproj)
wtk_albers = gpd.read_file("data/shapefiles/wtk_points_albers.shp")

# Rasterize both geographic and projected WTK points
src = "data/shapefiles/wtk_points.shp"
dst = "data/rasters/test.tif"
xmin, ymin, xmax, ymax = wtk.geometry.total_bounds
rasterize(src=src, dst=dst, attribute="gid", resolution=0.02, epsg=4326,
          extent=[xmin, ymin, xmax, ymax], overwrite=False)

src = "data/shapefiles/wtk_points_albers.shp"
dst = "data/rasters/test_albers.tif"
xmin, ymin, xmax, ymax = wtk_albers.geometry.total_bounds
rasterize(src=src, dst=dst, attribute="gid", resolution=, epsg=102008,
          extent=[xmin, ymin, xmax, ymax], overwrite=True)




# A sample shapefile layer
datapath = os.path.expanduser("~/Box/WETO 1.2/data")
samplepath = glob(os.path.join(datapath, "BLM", "fedlands_blm.shp"))[0]
sample = gpd.read_file(samplepath)


# 1) Create a GeoTiff from a shapefile.
shp =  "/Users/twillia2/Box/WETO 1.2/data/BLM/fedlands_blm.shp"
trgt = "data/fedlands.tif"

#



# 1) Create an HDF5 data set that lines up with the WTK CONUS data set




















# Let's clip by extent
xmin, ymin, xmax, ymax = wtk.geometry.total_bounds
test = sample.cx[xmin:xmax, ymin:ymax]

# We'll probably want to store these as both a geotiff and an hdf5 file
# GeoTiff




# Save to shapefile
test.to_file("data/test.shp")

# I can't infer the proper resolution, the projection created slightly
# different spacing...need to discuss this
src = "data/test.shp"
dst = "data/test_albers2.shp"

reproject_polygon(src, dst, tproj)


# WTK Points
wtk.to_file("data/wtk.shp")
src = "data/wtk.shp"
dst = "data/wtk_albers.shp"
reproject_point(src, dst, tproj)





# Is the project resolution consistent?
albers = gpd.read_file("data/test_albers2.shp")

src = "data/test_albers2.shp"
dst = "data/test_albers.tif"
rasterize(src=src, dst=dst, attribute="raster_val", resolution=0.02, epsg=4326,
          extent=[xmin, ymin, xmax, ymax], overwrite=False)


test.plot(column="state")
