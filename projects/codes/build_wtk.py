# -*- coding: utf-8 -*-
"""
Create WTK point/raster datas sets

Get the original points from the eagle datasets if possible.

Created on Fri Jan 10 15:01:16 2020

@author: twillia2
"""
import os
from geofeather import from_geofeather, to_geofeather
import geopandas as gpd
import pandas as pd
from functions import to_geo, rasterize, reproject_point

# I might share this on github
os.chdir(os.path.expanduser("~/github/weto"))

def build_wtk():

    # 1) Create Grids of WTK GIDS
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

    # Reproject to NA Alber Equal Area Conic
    tproj = ("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 " +
             "+x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 " +
             "+units=m +no_defs")
    if not os.path.exists("data/shapefiles/wtk_points_albers.shp"):
        reproject_point(src="data/shapefiles/wtk_points.shp",
                        dst="data/shapefiles/wtk_points_albers.shp",
                        tproj=tproj)
    wtk_albers = gpd.read_file("data/shapefiles/wtk_points_albers.shp")

    # Rasterize both geographic and projected WTK points - what resolution?
    src = "data/shapefiles/wtk_points.shp"
    dst = "data/rasters/test.tif"
    if not os.path.exists(dst):
        xmin, ymin, xmax, ymax = wtk.geometry.total_bounds
        rasterize(src=src, dst=dst, attribute="gid", resolution=0.02, epsg=4326,
                  extent=[xmin, ymin, xmax, ymax], overwrite=False)

    src = "data/shapefiles/wtk_points_albers.shp"
    dst = "data/rasters/test_albers.tif"
    if not os.path.exists(dst):
        xmin, ymin, xmax, ymax = wtk_albers.geometry.total_bounds
        rasterize(src=src, dst=dst, attribute="gid", resolution=49.5,
                  epsg=102008, extent=[xmin, ymin, xmax, ymax], overwrite=True)


if __name__ == "__main__":
    build_wtk()
