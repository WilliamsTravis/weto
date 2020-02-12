#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for WETO Non-Construction Costs input layer generation.

Created on Wed Jan  8 09:05:39 2020

@author: twillia2
"""
import dask.array as da
import os
import shutil
import geopandas as gpd
import h5py
import numpy as np
import rasterio
import subprocess as sp
import sys
import xarray as xr
from fiona.errors import DriverError
from multiprocessing import Pool
from osgeo import gdal, ogr, osr
from shapely.geometry import Point
from tqdm import tqdm

# Use GDAL exceptions
gdal.UseExceptions()

# VARIABLES
ALBERS_PROJ4 = ("+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 " +
                "+y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# Wind Toolkit Point Coordinates
# WTK = from_geofeather("data/shapefiles/wtk_points.feather")

@da.as_gufunc(signature="(i)->(i)", output_dtypes=int, vectorize=True,
              allow_rechunk=True)
def gmask(array, navalues=[-9999, 0]):
    """Create a mask of 1's and 0 out of an array using non values"""
    # So multiple values can be masked out
    if not isinstance(navalues, (list, tuple, np.ndarray)):
        navalues = [navalues]

    # Create mask
    array = array.copy()
    array[np.isin(array, navalues)] = np.nan
    array = array * 0 + 1
    array[np.isnan(array)] = 0

    return array


def to_xr(array, template):
    """Take a numpy or dask array and a georeferenced template xarray dataset
    to create a new georeferenced xarray dataset.

    Parameters
    ----------
    array : TYPE
        DESCRIPTION.
    template : TYPE
        DESCRIPTION.

    Returns
    -------
    dataset : TYPE
        DESCRIPTION.
    """

    # Get coordinates dims and attributes from the template
    dims = template.coords.dims
    coords = [template.coords[dim].data for dim in dims]
    attrs = template.attrs

    # Create a new data array and dataset
    darray= xr.DataArray(data=array, coords=coords, dims=dims)
    dataset = xr.Dataset(data_vars={'value': darray}, attrs=attrs)

    return dataset



# FUNCTIONS
def mask(array, navalues=-9999):
    """Create a mask of 1's and nans out of a numpy array"""

    # So that multiple values can be masked out
    if not isinstance(navalues, (list, tuple, np.ndarray)):
        navalues = [navalues]

    # Create mask
    mask = array.copy()
    mask[np.isin(mask, navalues)] = np.nan
    mask = mask * 0 + 1

    return mask


def rasterize(src, dst, attribute, epsg, transform, height, width,
              navalue=-9999, all_touch=False, dtype=gdal.GDT_Float32,
              overwrite=False):
    """
    Use GDAL Rasterize to rasterize a shapefile stored on disk and write
    outputs to a file.

    Parameters
    ----------
    src : str
        File path for the source file to rasterize.
    dst : str
        Destination path for the output raster.
    attribute : str
        Attribute name being rasterized.
    epsg : int
        EPSG Code associated with target coordinate reference system.
    transform : list | tuple | array
        Geometric affine transformation:
            (x-min, x-resolution, x-rotation, y-max, y-rotation, y-resoltution)
    height : int
        Number of y-axis grid cells.
    width : int
        Number of x-axis grid cells.
    na : int | float
        The value to assign to non-value grid cells. (defaults to -99999)
    all_touch : boolean
        Wether or not to associate vector values with all intersecting grid
        cells. (defaults to False)
    overwrite : boolean

    Returns
    -------
    None

    # Things to do:
        1) catch exceptions
        2) print stdout progress:
            https://gis.stackexchange.com/questions/237479/using-callback-with-
            python-gdal-rasterizelayer ?
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
    xmin, xres, xrot, ymax, yrot, yres = transform
    xs = [xmin + xres * i  for i in range(width)]
    ys = [ymax + yres * i  for i in range(height)]
    nx = len(xs)
    ny = len(ys)

    # Create the target raster layer
    driver = gdal.GetDriverByName("GTiff")
    trgt = driver.Create(dst, nx, ny, 1, dtype)
    trgt.SetGeoTransform((xmin, xres, xrot, ymax, yrot, yres))  # xres again?

    # Add crs
    refs = osr.SpatialReference()
    refs.ImportFromEPSG(epsg)
    trgt.SetProjection(refs.ExportToWkt())

    # Set no value
    band = trgt.GetRasterBand(1)
    band.SetNoDataValue(navalue)

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
    array[array == navalue] = np.nan

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
    src_file = None
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
    src_file = None
    trgt_file = None



def reproject_raster(src, dst, s_srs=None, t_srs=None, res=None, extent=None,
                     template=None, dtype="Float32", overwrite=False):
    """
    Warp a raster to a new geometry

    Parameters
    ----------
    src : str
        Path to source raster file.
    dst : str
        Path to target raster file.
    s_srs : str
        Source coordinate reference system. Requires a proj4 string for now.
    t_srs : str
        Target coordinate reference_system. Requires a proj4 string for now.
    res : int |float
        Target spatial resolution.
    extent : list-like
        Spatial extent of output raster in this order: (xmin, ymin, xmax, ymax)
    overwrite : boolean

    Returns
    -------
    None.

    sample arguments
        src = "/projects/rev/data/conus/_windready_conus.tif"
        dst = "/scratch/twillia2/weto/data/rasters/albers/exclusions.tif"
        template = "/scratch/twillia2/weto/data/rasters/albers/cost_codes.tif"
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


    # If a template is provided, use its geometry  for target figures
    if template:
        temp = rasterio.open(template)
        extent = list(temp.bounds)
        transform = temp.transform  # rasterio order is different that GDAL
        t_srs = temp.crs
        res = transform[0]

    # Get source srs
    source = rasterio.open(src)
    s_srs = source.crs

    # The data type needs to be a gdal object, but that's an odd thing to know
    dtypes ={
             "byte": gdal.GDT_Byte,
             "int16": gdal.GDT_Int16,
             "int32": gdal.GDT_Int32,
             "float32": gdal.GDT_Float32,
             "float64": gdal.GDT_Float64,
             "cint16": gdal.GDT_CInt16,
             "cint32": gdal.GDT_CInt32,
             "cfloat32": gdal.GDT_CFloat32,
             "cfloat64": gdal.GDT_CFloat64,
             }


    # Check Options: https://gdal.org/python/osgeo.gdal-module.html#WarpOptions
    try:
        ops = gdal.WarpOptions(format="GTiff",
                               outputBounds=extent,
                               xRes=res,
                               yRes=res,
                               srcSRS=s_srs,
                               dstSRS=t_srs,
                               resampleAlg="near",
                               outputType=dtype)
    except Exception as error:  # GDAL exceptions?
        print("Options not available: ")
        print(error)

    # Call
    ds = gdal.Warp(dst, src, options=ops)
    del ds


# def shp_to_h5(shp, hdf, attribute, mode="w"):
#     """Assign values from a shapefile of polygons to the WTK point dataset
#     return a geopandas geodataframe and write or append to an HDF5 file.

#     Parameters
#     ----------
#     shp : geopandas.geodataframe.GeoDataFrame | str
#         Either a GeoPandas GeoDataFrame or a path as a string to a shapefile
#     attribute : str
#         The attribute of the shapefile from which to assign values to the WTK
#         point data set.

#     Returns
#     -------
#     geopandas.geodataframe.GeoDataFrame

#     Mode Notes
#     ----------
#         r         : 	Readonly, file must exist
#         r+        : 	Read/write, file must exist
#         w         : 	Create file, truncate if exists
#         w- or x   : 	Create file, fail if exists
#         a         : 	Read/write if exists, create otherwise (default)

#     I want this to either create a new HDF5 file of cost layer values at each
#     wind toolkit point, or append to an existing hdf dataset with wtk points.
#     Grant will have the best idea of what is needed to incorporate our changes
#     into reV.

#     I should do multiple attributes at once.

#     Also, this drops non-present values...

#     It's a start
#     """

#     # Is this file in memory or on disk?
#     if isinstance(shp, str):
#         try:
#             shp = gpd.read_file(shp)
#         except DriverError:
#             print("Shapefile doesn't exist.")

#     # Extract just the requested attribute
#     shp = shp[["geometry", attribute]]
#     wtk2 = gpd.sjoin(WTK, shp)

#     # HDF5 doesn't like type 'O'
#     dtype = wtk2[attribute].dtype

#     if mode == "w":
#         with h5py.File(hdf, mode) as file:
#             if dtype is np.dtype('O'):
#                 file.create_dataset(name=attribute, data=wtk2[attribute],
#                                     dtype=h5py.string_dtype(encoding='utf-8'))
#             else:
#                 file.create_dataset(name=attribute, data=wtk2[attribute])
#             file.create_dataset(name="lat", data=wtk2["lat"])
#             file.create_dataset(name="lon", data=wtk2["lon"])
#     elif mode == "a":
#         with h5py.File(hdf, mode) as file:
#             if dtype is np.dtype('O'):
#                 file.create_dataset(name=attribute, data=wtk2[attribute],
#                                     dtype=h5py.string_dtype(encoding='utf-8'))
#             else:
#                 file.create_dataset(name=attribute, data=wtk2[attribute])
#     else:
#         print("Sorry, I haven't figured mode=" + mode + " out yet.")


# To split our big file into 100 smaller ones...
def split_extent(raster_file, n=100):
    """Split a raster files extent into n extent pieces"""

    # Get raster geometry
    rstr = rasterio.open(raster_file)
    geom = rstr.get_transform()
    ny = rstr.height
    nx = rstr.width
    xs = [geom[0] + geom[1] * i  for i in range(nx)]
    ys = [geom[3] + geom[-1] * i  for i in range(ny)]

    # Get number of chunks along each axis
    nc = np.sqrt(n)

    # Split coordinates into 10 pieces along both axes...
    xchunks = np.array(np.array_split(xs, nc))
    ychunks = np.array(np.array_split(ys, nc))

    # Get min/max of each coordinate chunk...
    sides = lambda x: [min(x), max(x)]
    xmap = map(sides, xchunks)
    ymap = map(sides, ychunks)
    xext = np.array([v for v in xmap])
    yext = np.array([v for v in ymap])

    # Combine these in this order [min, ymin, max, ymax]....
    extents = []
    for xex in xext:
        for yex in yext:
            extents.append([xex[0], yex[0], xex[1], yex[1]])

    return extents


def to_geo(data_frame, loncol="lon", latcol="lat", epsg=4326):
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
    data_frame["geometry"] = data_frame.apply(to_point, axis=1)
    gdf = gpd.GeoDataFrame(data_frame, geometry="geometry", crs=crs)

    return gdf


def to_raster(array, savepath, crs, geometry, navalue=-9999):
    """Takes in a numpy array and writes data to a GeoTiff.

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


def tileit(arg):
    """Use gdal to cut a raster into a smaller pieces"""

    # Separate arguments
    extent = arg[0]
    rfile = arg[1]
    chunk = arg[2]
    outfolder = arg[3]

    # Get everything in order
    extent = [str(e) for e in extent]
    chunk = "{:02d}".format(chunk)
    outbase =  os.path.basename(rfile).split(".")[0]
    outfile = os.path.join(outfolder, outbase + "_" + chunk + ".tif")

    # Let's not overwrite - check the file is good first actually
    if not os.path.exists(outfile):
        sp.call(["gdalwarp",
                 "-q",
                 "-te", extent[0], extent[1], extent[2], extent[3],
                 rfile,
                 outfile],
                stdout=sp.PIPE, stderr=sp.PIPE)

    return outfile


def tile_raster(raster_file, outfolder_path, ntiles, ncpu):
    """ Take a raster and write n tiles from it.

    Parameters
    ----------
    raster_file : str
        Path to a GeoTiff
    outfolder_path : str
        Path to a folder in which to store tiles. Will create if not present.
    ntiles : int
        Number of tiles to write.
    ncpu : int
        Number of cpus to use for processing.

    Returns
    -------
    None.


    sample args
    -----------
    raster_path = '/scratch/twillia2/weto/data/rasters/agcounty_product.tif'
    outfolder_path = '/scratch/twillia2/weto/data/rasters/chunks'
    ntiles=100
    ncpu=5
    """
    # Create tile folder
    base_name = os.path.splitext(raster_file)[0]
    out_folder = "_".join([base_name, "tiles"])
    os.makedirs(out_folder, exist_ok=True)

    # Get all of the extents needed to make n tiles
    extents = split_extent(raster_file, n=ntiles)

    # Wrap arguments into one object
    raster_files = np.repeat(raster_file, len(extents))
    chunknumbers = [i for i in range(len(extents))]
    out_folders = np.repeat(out_folder, len(extents))
    args = list(zip(extents, raster_files, chunknumbers, out_folders))

    # Run each
    with Pool(ncpu) as pool:
        tfiles = []
        for tfile in tqdm(pool.imap(tileit, args), total=len(extents),
                          position=0, file=sys.stdout):
            tfiles.append(tfile)

    return tfiles


def try_get(mapvals, val, errval=-9999):
    """Catching exceptions when vectorizing a dictionary mapping"""
    try:
        x = mapvals[val]
    except KeyError:
        x = errval
    return x


def vectorit(arg):
    """Map dictionary values from one raster to another"""
    # Get arguments
    infile = arg[0]
    outfile = arg[1]
    mapvals = arg[2]
    if not os.path.exists(outfile):
        ds = gdal.Open(infile)
        crs = ds.GetProjection()
        geom = ds.GetGeoTransform()
        array = ds.ReadAsArray()
        try:
            newarray = np.vectorize(try_get)(mapvals, array)
            to_raster(newarray, outfile, crs, geom, navalue=-9999)
        except Exception as e:
            print("\n")
            print(infile + ": ")
            print(e)
            print("\n")


def map_tiles(tilefiles, mapvals, outfolder, ncpu):
    """Take a list of tiled raster files, map value from a dictionary to
    a list of output raster files"""

    # Create the output paths
    os.makedirs(outfolder, exist_ok=True)
    outfiles = []
    for f in tilefiles:
        outfile = os.path.basename(f).replace("agcounty_product" ,"nlcd_index")
        outfiles.append(os.path.join(outfolder, outfile))

    # Bundle the arguments for vectorit (single function)
    dicts = [mapvals.copy() for i in range(len(outfiles))]
    args = list(zip(tilefiles, outfiles, dicts))

    # Run it
    with Pool(ncpu) as pool:
        for f in tqdm(pool.imap(vectorit, args), total=len(outfiles),
                      position=0, file=sys.stdout):
            pass

    # Return the output file paths
    return outfiles


def merge_tiles(tilefiles, outfile):
    """Take a list of paths to tiles raster files and merge them into a
    single raster.
    """
    gdal.merge


# CLASSES
class Data_Path:
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.expand_check()

    def join(self, file_path):
        return os.path.join(self.folder_path, file_path)

    def expand_check(self):
        if "~" in self.folder_path:
            self.folder_path =  os.path.expanduser(self.folder_path)
