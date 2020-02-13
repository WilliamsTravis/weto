"""
A set of functions for performing common spatial transformations using GDAL
bindings, geopandas, rasterio, and shapely.

Things to do:
    - Incorporate data types into these functions.
    - Continuously incorporate other GDAL functionality, too.
    - Some of these would be better placed in the spatial or utilities modules.
    - Create a file checking function and use it when writing new files. This
      could detect the file type, try to open it with the appropriate function,
      and raise an exception and delete it if it fails.
    - use **kwargs to include all available options. The check will still flag
      non-extant options. This will also make the rest of the gdal functions
      easier to implement.
    - add in creation options, as these are separate from function options.
      This is trickier though, because there can be several. Make this a list
      in an explicit "co" argument?
"""
import dask.array as da
import geopandas as gpd
import numpy as np
import os
import rasterio
import shutil
import subprocess as sp
import sys
from multiprocessing import Pool
from osgeo import gdal, ogr, osr
from shapely.geometry import Point
from tqdm import tqdm
import weto.dask_raster as dr


# CONSTANTS
GDAL_TYPES = {"GDT_Byte": "Eight bit unsigned integer",
              "GDT_CFloat32": "Complex Float32",
              "GDT_CFloat64": "Complex Float64",
              "GDT_CInt16": "Complex Int16",
              "GDT_CInt32": "Complex Int32",
              "GDT_Float32": "Thirty two bit floating point",
              "GDT_Float64": "Sixty four bit floating point",
              "GDT_Int16": "Sixteen bit signed integer",
              "GDT_Int32": "Thirty two bit signed integer",
              "GDT_UInt16": "Sixteen bit unsigned integer",
              "GDT_UInt32": "Thirty two bit unsigned integer",
              "GDT_Unknown": "Unknown or unspecified type"}

GDAL_TYPEMAP = {"byte": gdal.GDT_Byte,
                "cfloat32": gdal.GDT_CFloat32,
                "cfloat64": gdal.GDT_CFloat64,
                "cint16": gdal.GDT_CInt16,
                "cint32": gdal.GDT_CInt32,
                "float32": gdal.GDT_Float32,
                "float64": gdal.GDT_Float64,
                "int16": gdal.GDT_Int16,
                "int32": gdal.GDT_Int32,
                "uint16": gdal.GDT_UInt16,
                "uint32": gdal.GDT_UInt32,
                "unknown": gdal.GDT_Unknown}


# FUNCTIONS
def gdal_progress(complete, message, unknown):
    """A progress callback that recreates the gdal printouts."""

    # We don't need the message or unknown objects
    del message, unknown

    # Complete is a number between 0 and 1
    percent = int(complete * 100)

    # Between numeric printouts we need three dots
    dots = [[str(i) + d for d in ["2", "5", "8"]] for i in range(10)]
    dots = [int(l) for sl in dots for l in sl]

    # If divisible by ten, print the number
    if percent % 10 == 0 and percent != 0:
        print("{}".format(percent), end="")

    # If one of three numbers between multiples of 10, print a dot
    elif percent in dots:
        print(".", end="")

    return 1


def check_raster(file):
    """Check that an output raster opens, delete if it fails

    Parameters
    ----------
    file : str
        Path to raster file.

    Returns
    -------
    None.
    """

    try:
        file = gdal.Open(file)
    except:  # <--------------------------------------------------------------- Catch exceptions as they come up.
        print(file + " failed to open.")


def rasterize(src, dst, attribute, epsg, transform, height, width,
              navalue=-9999, all_touch=False, dtype=gdal.GDT_Float32,
              template=None, overwrite=False):
    """
    Use GDAL RasterizeLayer to rasterize a shapefile stored on disk and write
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
    dtype : str | gdal object
        GDAL data type. Can be a string or a gdal type object (e.g.
        gdal.GDT_Float32, "GDT_Float32", "float32"). Available GDAL data types
        and descriptions can be found in the GDAL_TYPES dictionary.
    overwrite : boolean

    Returns
    -------
    None

    # Things to do:
        1) Catch exceptions
        2) Progress callback
        3) Use a template as option.
        4) Use more than just EPSG (doesn't always work, also accept proj4)
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

    # if template:
        # Include as much as possible from the template

    # Use transform to derive coordinates and dimensions
    xmin, xres, xrot, ymax, yrot, yres = transform
    xs = [xmin + xres * i for i in range(width)]
    ys = [ymax + yres * i for i in range(height)]
    nx = len(xs)
    ny = len(ys)

    # Specifying data types shouldn't be so difficult
    if isinstance(dtype, str):
        dtype = dtype.lower().replace("gdt_", "")
        try:
            dtype = GDAL_TYPEMAP[dtype]
        except KeyError:
            print("\n'" + dtype + "' is not an available data type. "
                  "Choose a value from this list:")
            print(str(list(GDAL_TYPEMAP.keys())))

    # Create the target raster layer
    driver = gdal.GetDriverByName("GTiff")
    trgt = driver.Create(dst, nx, ny, 1, dtype)
    trgt.SetGeoTransform((xmin, xres, xrot, ymax, yrot, yres))

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

    # Alternately
    # try:
          # ops = gdal.RasterizeOptions(options=None,
          #                             format=None,
          #                             outputType=0,
          #                             creationOptions=None,
          #                             noData=None,
          #                             initValues=None,
          #                             outputBounds=None,
          #                             outputSRS=None,
          #                             transformerOptions=None,
          #                             width=None,
          #                             height=None,
          #                             xRes=None,
          #                             yRes=None,
          #                             targetAlignedPixels=False,
          #                             bands=None,
          #                             inverse=False,
          #                             allTouched=False,
          #                             burnValues=None,
          #                             attribute=None,
          #                             useZ=False,
          #                             layers=None,
          #                             SQLStatement=None,
          #                             SQLDialect=None,
          #                             where=None,
          #                             optim=None,
          #                             callback=None,
          #                             callback_data=None)
    # except Exception as error:  # GDAL exceptions?
    #     print("Options not available: ")
    #     print(error)



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
             raster values : numpy.ndarray
             affine transformation : tuple
                 (top left x coordinate, x resolution, row rotation,
                  top left y coordinate, column rotation, y resolution)),
            coordinate reference system : str
                 Well-Known Text format
    """

    # Open raster file and read in parts necessary for rewriting
    raster = gdal.Open(rasterpath)
    geometry = raster.GetGeoTransform()
    arrayref = raster.GetProjection()
    array = np.array(raster.GetRasterBand(band).ReadAsArray())
    raster = None

    # This helped for some old use-case, but might not be necessary
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
    This only handles ESRI Shapefiles at the moment, but can be written to
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
    except Exception:
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
    This only handles ESRI Shapefiles at the moment, but can be written to
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
    except Exception:
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


def split_extent(raster_file, n=100):
    """Split a raster files extent into n extent pieces."""

    # Get raster geometry
    rstr = rasterio.open(raster_file)
    geom = rstr.get_transform()
    ny = rstr.height
    nx = rstr.width
    xs = [geom[0] + geom[1] * i for i in range(nx)]
    ys = [geom[3] + geom[-1] * i for i in range(ny)]

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

    # Combine these in this order [xmin, ymin, xmax, ymax]....
    extents = []
    for xex in xext:
        for yex in yext:
            extents.append([xex[0], yex[0], xex[1], yex[1]])

    return extents


def tile_raster(raster_file, out_folder, ntiles, ncpu):
    """ Take a raster and write n tiles from it.

    Parameters
    ----------
    raster_file : str
        Path to a GeoTiff
    out_folder : str
        Path to a folder in which to store tiles. Will create if not present.
    ntiles : int
        Number of tiles to write.
    ncpu : int
        Number of cpus to use for processing.

    Returns
    -------
    None.
    """

    # Create the output folder
    if not out_folder:
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
        for tfile in tqdm(pool.imap(tile_single, args), total=len(extents),
                          position=0, file=sys.stdout):
            tfiles.append(tfile)

    return tfiles


def tile_single(arg):
    """Use gdal to cut a raster into a smaller pieces.

    Note:
        This is made for tile_raster and is not intuitive as a standalone.
        Add in a check to make sure each output file is good. Moving to
        a class method soon.
    """

    # Separate arguments
    extent = arg[0]
    rfile = arg[1]
    chunk = arg[2]
    outfolder = arg[3]

    # Get everything in order
    extent = [str(e) for e in extent]
    chunk = "{:02d}".format(chunk)
    outbase = os.path.basename(rfile).split(".")[0]
    outfile = os.path.join(outfolder, outbase + "_" + chunk + ".tif")

    # Let's not overwrite - use the warp function from below here instead
    if not os.path.exists(outfile):
        sp.call(["gdalwarp",
                 "-q",
                 "-te", extent[0], extent[1], extent[2], extent[3],
                 rfile,
                 outfile],
                stdout=sp.PIPE, stderr=sp.PIPE)

    return outfile


def to_geo(data_frame, loncol="lon", latcol="lat", epsg=4326):
    """Convert a Pandas DataFrame object to a GeoPandas GeoDataFrame object.

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


def to_raster(array, savepath, crs=None, geometry=None, template=None,
              dtype=gdal.GDT_Float32, tiled=False, compress=None,
              nadata=-9999):
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
    template : str
        Path to a raster file with desired target raster geometry, crs, and na
        value. This will overwrite other arguments provided for these
        parameters.
    dtype : str | gdal object
        GDAL data type. Can be a string or a gdal type object (e.g.
        gdal.GDT_Float32, "GDT_Float32", "float32"). Available GDAL data types
        and descriptions can be found in the GDAL_TYPES dictionary.
    compress : str
        A compression technique. Available options are "DEFLATE", "JPEG",
        "LZW"
    navalue : int | float
        The number used for non-values in the raster data set. Defaults to
        -9999.


    # To Do:
        - Write multi band rasters
        - Write dask arrays directly to disk
        - Use different drivers
        -
    """

    # Retrieve needed raster elements
    xpixels = array.shape[1]
    ypixels = array.shape[0]

    # Use a template file to extract affine transformation, crs, and na value
    if template:
        template_file = gdal.Open(template)
        geometry = template_file.GetGeoTransform()
        crs = template_file.GetProjection()
        band1 = template_file.GetRasterBand(1)
        nadata = band1.GetNoDataValue()

    # Get options here - not built out yet
    if compress:
        cops = ["compress=" + compress]

    # if this is a dask array...
    if isinstance(array, da.core.Array):  # <---------------------------------- Just use gdal here to be consistent.
        print("DASK ARRAY!")
        template_file = rasterio.open(template)
        transform = template_file.transform
        with dr.RasterioDataset(savepath, 'w',
                                driver='GTiff',
                                width=xpixels,
                                height=ypixels,
                                count=1,
                                crs=crs,
                                transform=transform,
                                dtype=dtype,
                                tiled=tiled,
                                nadata=nadata,
                                compress=compress) as file:
            da.store(array, file, lock=True)
    else:
        # Specifying data types shouldn't be so difficult
        if isinstance(dtype, str):
            dtype = dtype.lower().replace("gdt_", "")
            try:
                dtype = GDAL_TYPEMAP[dtype]
            except KeyError:
                print("\n'" + dtype + "' is not an available data type. "
                      "Choose a value from this list:")
                print(str(list(GDAL_TYPEMAP.keys())))

        # Create file
        driver = gdal.GetDriverByName("GTiff")
        image = driver.Create(savepath, xpixels, ypixels, 1, dtype,
                              options=cops)

        # Write raster data and attributes to file
        image.SetGeoTransform(geometry)
        image.SetProjection(crs)
        image.GetRasterBand(1).WriteArray(array)
        image.GetRasterBand(1).SetNoDataValue(nadata)


def warp(src, dst, s_srs=None, t_srs=None, res=None, extent=None,
         dtype="Float32", src_nodata=-9999., dst_nodata=-9999., template=None,
         overwrite=False, compress=None):
    """
    Warp a raster to a new geometry.

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
    dtype : str | gdal object
        GDAL data type. Can be a string or a gdal type object (e.g.
        gdal.GDT_Float32, "GDT_Float32", "float32"). Available GDAL data types
        and descriptions can be found in the GDAL_TYPES dictionary.
    src_nodata : int | float
        Numeric value in the input raster to interpret as non-data.
    dst_nodata : int | float
        Numeric value to assign as non-data in the output raster.
    template : str
        Path to a raster file with desired target raster geometry, crs,
        resolution, and extent values. This will overwrite other arguments
        provided for these parameters.
    overwrite : boolean
    compress : str
        A compression technique. Available options are "DEFLATE", "JPEG",
        "LZW"

    Returns
    -------
    None.
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

    # Specifying data types shouldn't be so difficult
    if isinstance(dtype, str):
        dtype = dtype.lower().replace("gdt_", "")
        try:
            dtype = GDAL_TYPEMAP[dtype]
        except KeyError:
            print("\n'" + dtype + "' is not an available data type. "
                  "Choose a value from this list:")
            print(str(list(GDAL_TYPEMAP.keys())))

    # Compress
    if compress:
        co = ["COMPRESS=" + compress]
    else:
        co = None

    # Check Options: https://gdal.org/python/osgeo.gdal-module.html#WarpOptions
    try:
        ops = gdal.WarpOptions(format="GTiff",
                               outputBounds=extent,
                               xRes=res,
                               yRes=res,
                               srcSRS=s_srs,
                               dstSRS=t_srs,
                               resampleAlg="near",
                               outputType=dtype,
                               srcNodata=src_nodata,
                               dstNodata=dst_nodata,
                               callback=gdal_progress,
                               creationOptions=co)

    except Exception as error:  # GDAL exceptions?
        print("Options not available: ")
        print(error)

    # Call
    print("Processing " + dst + " :")
    ds = gdal.Warp(dst, src, options=ops)
    del ds


# CLASSES
class Data_Path:
    """Data_Path joins a root directory path to data file paths."""

    def __init__(self, folder_path):
        """Initialize Data_Path."""

        self.folder_path = folder_path
        self._expand_check()

    def join(self, file_path):
        """Join a file path to the root directory path"""

        return os.path.join(self.folder_path, file_path)

    def _expand_check(self):

        # Expand the user path if a tilda is present in the root folder path.
        if "~" in self.folder_path:
            self.folder_path = os.path.expanduser(self.folder_path)


class Map_Values:
    """Map a set of keys from an input raster (or rasters) to values in an
    output raster (or rasters) using a dictionary of key-value pairs."""

    def __init__(self, val_dict, err_val=-9999):
        """Initialize Map_Values.

        Parameters
        ----------
        val_dict : dict
            A dictionary of key-value pairs
        errval : int | float
            A value to assign where there are no matching keys in val_dict.
        """
        self.val_dict = val_dict
        self.err_val = err_val

    def map_file(self, src, dst):
        """Take an input raster file, map values from a dictionary to an output
        raster file.

        Parameters
        ----------
        src : str
            Path to the input raster file.
        dst : str
            Path to the output raster file. Directory will be created if it
            does not exist.

        Returns
        -------
        None.
        """

        # Create the output path
        out_folder = os.path.dirname(dst)
        os.makedirs(out_folder, exist_ok=True)

        # Bundle the arguments for map_single (single function)
        args = list(src, dst, self.val_dict)

        # Run it
        self._map_single(args)

    def map_files(self, src_files, out_folder, ncpu):
        """Take a list of tiled raster files, map values from a dictionary to
        a list of output raster files.

        Parameters
        ----------
        src_files : list-like
            A list of paths to raster files.
        outfolder : str
            A path to a target directory to store output files. Will be
            created if it does not exist.
        ncpu : int
            The number of cpus to use for multiprocessing.

        Returns
        -------
        outfiles : list
            A list of paths to output files.
        """

        # Create the output paths
        os.makedirs(out_folder, exist_ok=True)
        dst_files = []
        for file in src_files:
            dst_file = os.path.basename(file)
            dst_files.append(os.path.join(out_folder, dst_file))

        # Bundle the arguments for map_single (single function)
        dicts = [self.val_dict.copy() for i in range(len(dst_files))]
        args = list(zip(src_files, dst_files, dicts))

        # Run it
        with Pool(ncpu) as pool:
            for _ in tqdm(pool.imap(self._map_single, args), position=0,
                          total=len(dst_files), file=sys.stdout):
                pass

        # Return the output file paths
        return dst_files

    def _map_single(self, arg):
        """Map dictionary values from one raster file to another.

        Parameters
        ----------
        arg : list-like
            A list containing an input raster file path, and output raster file
            path and a dictionary (bundled for multiprocessing).

        Returns
        -------
        None.
        """

        # Get arguments
        src = arg[0]
        dst = arg[1]
        val_dict = arg[2]

        # Try to map values from the mapvals dictionary to a new raster
        if not os.path.exists(dst):
            ds = gdal.Open(src)
            crs = ds.GetProjection()
            geom = ds.GetGeoTransform()
            array = ds.ReadAsArray()
            try:
                new_array = np.vectorize(self._map_try)(val_dict, array)
                to_raster(new_array, dst, crs, geom, navalue=-9999)
            except Exception as error:
                print("\n")
                print(src + ": ")
                print(error)
                print("\n")

    def _map_try(self, val_dict, key):
        """Use a key to return a dictionary value, return a specified value for
        exceptions.

        Parameters
        ----------
        val_dict : dict
            A dictionary of values.
        key : str | int | float
            A key that corresponds to a value in val_dict.

        Returns
        -------
        x : int | float
            The value from val_dict corresponding with the key.
        """

        # Try to retrieve a value with the key
        try:
            x = val_dict[key]

        # Return a specified error value if the key is not present
        except KeyError:
            x = self.err_val

        return x


# SCRATCH
def gdal_options(module="warp", **kwargs):
    """Capture any availabe option for gdalwarp."""

    # All available options
    modules = [m for m in gdal.__dict__ if "Options" in m and "_" not in m]

    # Get the module option method associated with module
    # ...

    # Return the appropriate options object
    try:
        ops = gdal.WarpOptions(**kwargs)
        return ops
    except TypeError:  # GDAL exceptions?
        print("One or more of the Warp options provided are not available: ")
        docs = "\n".join(gdal.WarpOptions.__doc__.split("\n")[1:])
        print(docs)
