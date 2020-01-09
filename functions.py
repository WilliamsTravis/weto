#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for WETO Non-Construction Costs input layer generation.

Created on Wed Jan  8 09:05:39 2020

@author: twillia2
"""
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import shutil
from glob import glob
from osgeo import gdal, ogr, osr
from shapely.geometry import Point


# FUNCTIONS
def to_geo(df, loncol="lon", latcol="lat", epsg=4326):
    """ Convert a Pandas DataFrame object to a GeoPandas GeoDataFrame object.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        A pandas data frame with latitude and longitude coordinates.
    loncol : str
        The name of the longitude column.
    latcol : str
        The name of the latitude column.
    epsg : int
        EPSG code associated with the Coordinate Reference System.

    Returns
    -------
    geopandas.geodataframe.GeoDataFrame
        A GeoPandas GeoDataFrame object.
    """
    crs = {"init": "epsg:{}".format(epsg)}
    to_point = lambda x: Point((x[loncol], x[latcol]))
    df["geometry"] = df.apply(to_point, axis=1)
    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=crs)

    return gdf


def rasterize(src, dst, attribute, resolution, epsg, extent, na=-9999,
              all_touch=False, overwrite=False):
    """Use GDAL Rasterize to rasterize a shapefile stored on disk and write
    outputs to a file.

    Parameters
    ----------
    src : str
        File path for the source file to rasterize.
    dst : str
        Destination path for the output raster.
    attribute : str
        Attribute name being rasterized.
    resolution : int | float
        Desired grid cell resolution.
    epsg : int
        EPSG code associated with the data's coordinate reference system.
    extent : list | list-like object
        Target geographic extent of output raster in this order:
        [xmin, ymin, xmax, ymax]
    na : int | float
        The value to assign to non-value grid cells. (defaults to -99999)
    all_touch : boolean
        Wether or not to associate vector values with all intersecting grid
        cells. (defaults to False)
    overwrite : boolean
    """
    # Overwrite existing file
    if os.path.exists(dst):
        if overwrite:
            if os.path.isfile(dst):
                os.remove(dst)
            else:
                shutil.rmtree(dst)
        else:
            print(dst + " exists, use overwrite=True to replace this file.")
            return

    # Open shapefile, retrieve the layer
    src_data = ogr.Open(src)
    layer = src_data.GetLayer()

    # Use transform to derive coordinates and dimensions
    xmin = extent[0]
    ymin = extent[1]
    xmax = extent[2]
    ymax = extent[3]

    # Create the target raster layer
    cols = int((xmax - xmin)/resolution)
    rows = int((ymax - ymin)/resolution) + 1
    trgt = gdal.GetDriverByName("GTiff").Create(dst, cols, rows, 1,
                                gdal.GDT_Float32)
    trgt.SetGeoTransform((xmin, resolution, 0, ymax, 0, -resolution))

    # Add crs
    refs = osr.SpatialReference()
    refs.ImportFromEPSG(epsg)
    trgt.SetProjection(refs.ExportToWkt())

    # Set no value
    band = trgt.GetRasterBand(1)
    band.SetNoDataValue(na)

    # Set options
    if all_touch is True:
        ops = ["-at", "ATTRIBUTE=" + attribute]
    else:
        ops = ["ATTRIBUTE=" + attribute]

    # Finally rasterize
    gdal.RasterizeLayer(trgt, [1], layer, options=ops)

    # Close target an source rasters
    del trgt
    del src_data


def read_raster(rasterpath, band=1, navalue=-9999):
    """Converts a raster file on disk into a numpy array along with
    spatial features needed to write results to a raster file.

    Parameters
    ----------
    rasterpath : str
        Path to a raster file.
    band : int
        The band number desired.
    navalue : int | float
        The number used for non-values in the raster data set

    Returns
    -------
        tuple:
             raster values (numpy.ndarray),
             affine transformation (tuple: (top left x coordinate,
                                            x resolution, row rotation,
                                            top left y coordinate,
                                            column rotation, y resolution)
                                   ),
            coordinate reference system (str: Well-Known Text format)
    """
    raster = gdal.Open(rasterpath)
    geometry = raster.GetGeoTransform()
    arrayref = raster.GetProjection()
    array = np.array(raster.GetRasterBand(band).ReadAsArray())
    raster = None
    array = array.astype(float)
    if np.nanmin(array) < navalue:
        navalue = np.nanmin(array)
    array[array==navalue] = np.nan

    return(array, geometry, arrayref)


def reproject_polygon(src, dst, tproj):
    """Reproject a shapefile of polygons and write results to disk. Recreates
    this GDAL command:
        ogr2ogr -s_srs <source_projection> -t_srs <target_projection> dst src

    Parameters
    ----------
    src : str
        Path to source shapefile.
    dst : str
        Path to target file.
    tproj (int | str):
        Target coordinate projection system as an epsg code or proj4 string.
        Sometimes EPSG codes aren't available to GDAL installations, but
        they're easier to use when they are so this will try both.

    Note
    ----
    This only handles ESRI Shapefiles at the moment, but will be written to
    handle any available driver.
    """

    # Create Shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")

    # Source reference information
    src_file = driver.Open(src)
    src_layer = src_file.GetLayer()
    src_srs = src_layer.GetSpatialRef()
    src_defn = src_layer.GetLayerDefn()

    # Target reference information
    trgt_srs = osr.SpatialReference()
    try:
        trgt_srs.ImportFromEPSG(tproj)
    except:
        trgt_srs.ImportFromProj4(tproj)

    # The transformation equation
    transform = osr.CoordinateTransformation(src_srs, trgt_srs)

    # Target file and layer
    if os.path.exists(dst):
        driver.DeleteDataSource(dst)
    trgt_file = driver.CreateDataSource(dst)
    trgt_layer = trgt_file.CreateLayer('', trgt_srs, ogr.wkbMultiPolygon)

    # Add Fields
    for i in range(0, src_defn.GetFieldCount()):
        defn = src_defn.GetFieldDefn(i)
        trgt_layer.CreateField(defn)

    # Get the target layer definition
    trgt_defn = trgt_layer.GetLayerDefn()

    # You have to reproject each feature
    src_feature = src_layer.GetNextFeature()
    while src_feature:
        # Get geometry
        geom = src_feature.GetGeometryRef()

        # Reproject geometry
        geom.Transform(transform)

        # Create target feature
        trgt_feature = ogr.Feature(trgt_defn)
        trgt_feature.SetGeometry(geom)
        for i in range(0, trgt_defn.GetFieldCount()):
            trgt_feature.SetField(trgt_defn.GetFieldDefn(i).GetNameRef(),
                                  src_feature.GetField(i))

        # Add feature to target file
        trgt_layer.CreateFeature(trgt_feature)

        # Close current feature
        trgt_feature = None

        # Get the next feature
        src_feature = src_layer.GetNextFeature()

    # Close both shapefiles
    src_file= None
    trgt_file = None


def reproject_point(src, dst, tproj):
    """Reproject a shapefile of points and write results to disk. Recreates
    this GDAL command:

        ogr2ogr -s_srs <source_projection> -t_srs <target_projection> dst src

    Parameters
    ----------
    src : str
        Path to source shapefile.
    dst : str
        Path to target file.
    tproj (int | str):
        Target coordinate projection system as an epsg code or proj4 string.
        Sometimes EPSG codes aren't available to GDAL installations, but
        they're easier to use when they are so this will try both.

    Note
    ----
    This only handles ESRI Shapefiles at the moment, but will be written to
    handle any available driver.
    """

    # Create Shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")

    # Source reference information
    src_file = driver.Open(src)
    src_layer = src_file.GetLayer()
    src_srs = src_layer.GetSpatialRef()
    src_defn = src_layer.GetLayerDefn()

    # Target reference information
    trgt_srs = osr.SpatialReference()
    try:
        trgt_srs.ImportFromEPSG(tproj)
    except:
        trgt_srs.ImportFromProj4(tproj)

    # The transformation equation
    transform = osr.CoordinateTransformation(src_srs, trgt_srs)

    # Target file and layer
    if os.path.exists(dst):
        driver.DeleteDataSource(dst)
    trgt_file = driver.CreateDataSource(dst)
    trgt_layer = trgt_file.CreateLayer('', trgt_srs, ogr.wkbPoint)

    # Add Fields
    for i in range(0, src_defn.GetFieldCount()):
        defn = src_defn.GetFieldDefn(i)
        trgt_layer.CreateField(defn)

    # Get the target layer definition
    trgt_defn = trgt_layer.GetLayerDefn()

    # You have to reproject each feature
    src_feature = src_layer.GetNextFeature()
    while src_feature:
        # Get geometry
        geom = src_feature.GetGeometryRef()

        # Reproject geometry
        geom.Transform(transform)

        # Create target feature
        trgt_feature = ogr.Feature(trgt_defn)
        trgt_feature.SetGeometry(geom)
        for i in range(0, trgt_defn.GetFieldCount()):
            trgt_feature.SetField(trgt_defn.GetFieldDefn(i).GetNameRef(),
                                  src_feature.GetField(i))

        # Add feature to target file
        trgt_layer.CreateFeature(trgt_feature)

        # Close current feature
        trgt_feature = None

        # Get the next feature
        src_feature = src_layer.GetNextFeature()

    # Close both shapefiles
    src_file= None
    trgt_file = None


def to_raster(array, savepath, crs, geometry, navalue=-9999):
    """Takes in a numpy.ndarray and writes data to a GeoTiff.

    Parameters
    ----------
    array : numpy.ndarray
        Numpy array to write to raster file.
    savepath : str
        Path to the target raster file.
    crs : str
        Coordinate reference system in Well-Known Text format.
    geometry : tuple
        Affine transformation information in this order:
            (top left x coordinate, x resolution, row rotation,
            top left y coordinate, column rotation, y resolution)
    navalue : int | float
        The number used for non-values in the raster data set.
    """
    # Retrieve needed raster elements
    geometry = geometry
    crs = crs
    xpixels = array.shape[1]
    ypixels = array.shape[0]

    # This helps sometimes
    savepath = savepath.encode('utf-8')

    # Create file
    driver = gdal.GetDriverByName("GTiff")
    image = driver.Create(savepath, xpixels, ypixels, 1, gdal.GDT_Float32)

    # Save raster and attributes to file
    image.SetGeoTransform(geometry)
    image.SetProjection(crs)
    image.GetRasterBand(1).WriteArray(array)
    image.GetRasterBand(1).SetNoDataValue(navalue)
