# -*- coding: utf-8 -*-
"""
Functions to help determine soil categories and depth to rock bottom of CONUS using the
gridded SSURGO data set.


Some manual steps are currently needed:
    1) download ESRI GeoDatabase from:
        https://nrcs.app.box.com/v/soils/folder/94124173798 (conus) or
        https://nrcs.app.box.com/v/soils/folder/94128402340 (states)
    2) Extract the MapunitRaster from the geodatabase using:
        ESRI or
        https://github.com/r-barnes/ArcRasterRescue (rename to .tif)
    3) If you went for the rescuer route, the edges are a bit off:
       change all values less than 0 to the navalue, and write it back.
        
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016MS000686

Created on Tue Mar 17 12:57:29 2020

@author: travis
"""

import os
import subprocess as sp

import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio 
import xarray as xr

from dask.distributed import Client
from gdalmethods import Map_Values
from pyproj import Proj, transform
from urllib.request import urlretrieve


GSSURGO_URLS = {
    "conus": "https://nrcs.app.box.com/v/soils/folder/94124173798",
    "state": "https://nrcs.app.box.com/v/soils/folder/94128402340"
    }

GNATSO_URLS = {
    "conus": "https://nrcs.app.box.com/v/soils/folder/83297297479",
    "state": "https://nrcs.app.box.com/v/soils/folder/84152602915"
    }

RESCUE_INSTRUCTIONS = open(
    os.path.expanduser(
        "~/github/weto/projects/soil/raster_rescue_instructions.txt")
    ,"r").read()

def fix_mukey(mukey_tif):
    """The mukey tiff might come out with some odd values around the edges,
    which makes viewing it a bit tricky. This fixes that."""

    mukey_tif = os.path.expanduser(mukey_tif)

    # For the profile
    profile = rasterio.open(mukey_tif).profile

    # For the data array
    mukey = xr.open_rasterio(mukey_tif, chunks=(1, 5000, 5000))

    # Change everything under 0 to the nan value
    array = mukey.data
    array[array < 0] = profile["nodata"]
    with Client():
        array = array.compute()

    # save
    with rasterio.Env():
        with rasterio.open(mukey_tif, "w", **profile) as file:
            file.write(array[0].astype(rasterio.int32), 1)


def get_gssurgo(state=None, dst="./"):
    """
    Download and translate the Gridded Soil Survey Geographic Data Set from
    the National Resource Conservation Service. If a state acronym is provided
    this will download only data for that state.

    Parameters
    ----------
    state : str, optional
        Acronym of a US state. The default is None.
    dst : str, optional
        Path to target directly for storing gSSURGO. The default is "./".

    Returns
    -------
    None.
    """

    # Build URL
    if state:
        state_name = state.upper()
        base_url = GSSURGO_URLS["state"]
        url = os.path.join(base_url, "gSSURGO_" + state_name + ".gdb.zip")
    else:
        base_url = GSSURGO_URLS["conus"]
        url = os.path.join(base_url, "gSSURGO.gdb.zip")

    # Download file
    filename = os.path.join(dst, os.path.basename(url))
    try:
        urlretrieve(url, filename)
    except:
        print("Getting gSSURO. I actually haven't written this in yet. It "
              "seems to require a Box SDK API, but the authentication process "
              "is stupid since this is public. \n For states go to " +
              GSSURGO_URLS["state"] + ". \n For CONUS go to " +
              GSSURGO_URLS["conus"] + ".\n")


def get_gnatsgo(state=None, dst="./"):
    """
    Download and translate the Gridded Soil Survey Geographic Data Set from
    the National Resource Conservation Service. If a state acronym is provided
    this will download only data for that state.

    Parameters
    ----------
    state : str, optional
        Acronym of a US state. The default is None.
    dst : str, optional
        Path to target directly for storing gSSURGO. The default is "./".

    Returns
    -------
    None.
    """

    # Build URL
    if state:
        state_name = state.upper()
        base_url = GNATSO_URLS["state"]
        url = os.path.join(base_url, "gSSURGO_" + state_name + ".gdb.zip")
    else:
        base_url = GNATSO_URLS["conus"]
        url = os.path.join(base_url, "gSSURGO.gdb.zip")

    # Download file
    filename = os.path.join(dst, os.path.basename(url))
    try:
        urlretrieve(url, filename)
    except:
        print("Getting gNATSGO. I actually haven't written this in yet. It "
              "seems to require a Box SDK API, but the authentication process "
              "is stupid since this is public. \n For states go to " +
              GNATSO_URLS["state"] + ". \n For CONUS go to " +
              GNATSO_URLS["conus"] + ".\n")

def map_variable(gdb_path, mukey_path, variable, dst):
    """
    Create a dataset of a SSURGO soil variable at all available horizons. The
    dataset will contain n layers corresponding to the maximum number of
    available soil horizons across all cells.

    Parameters
    ----------
    gdb_path : str
        Path to ESRI gSSRUGO Geodatabase
    mukey_path : str
        Path to a raster of the gSSURGO's 10 m map unit key raster.
    variable : str
        Name of the soil variable to be mapped.
    dst : str
        Path to destination file.

    Notes
    -----
    Variable descriptions:
        https://data.nal.usda.gov/system/files/SSURGO_Metadata_-_Table_Column_Descriptions.pdf#page=81

    Units:
        https://jneme910.github.io/CART/chapters/Soil_Propert_List_and_Definition

    Sample Arguments
    ----------------
    mukey_path = "~/data/weto/soil/mukey_de.tif"
    gdb_path = "~/data/weto/soil/gNATSGO_DE.gdb"
    dst = "~/data/weto/soil/brockdepmin_de.tif"
    variable = "brockdepmin"
    """

    # Expand user path
    mukey_path = os.path.expanduser(mukey_path)
    gdb_path = os.path.expanduser(gdb_path)
    dst = os.path.expanduser(dst)

    # Get the Map Unit Aggregated Attribute Table
    mukey = xr.open_rasterio(mukey_path, chunks=(1, 5000, 5000))
    umukeys = np.unique(mukey[:]).astype(str)
    muaggatt = gpd.read_file(gdb_path, layer="muaggatt")
    chorizon = gpd.read_file(gdb_path, layer="chorizon")
    cogeomordesc = gpd.read_file(gdb_path, layer="cogeomordesc")
    components = gpd.read_file(gdb_path, layer="component")
    components = pd.merge(chorizon, components, on="cokey")
    components = pd.merge(components, muaggatt, on="mukey")

    # Put the keys in front
    keys = [c for c in components.columns if "key" in c]
    others = [c for c in components.columns if "key" not in c]
    new_order = keys + others
    components = components[new_order]

    # Get the Horizon Table
    munits = muaggatt[["mukey", "muname"]]
    variable_df = pd.merge(variable_df, munits, on="mukey")
    variable_df = components[["mukey", "chkey", "muname", "hzname", "geomdesc",
                              "desgnmaster", "hzdept_r", "hzdepb_r", "hzthk_r",
                              "sandtotal_r", "silttotal_r", "claytotal_r", 
                              variable]]

    # Now, whats the best way to map these values
    val_dict = dict(zip(variable_df["mukey"].astype(int),
                        variable_df[variable]))
    mv = Map_Values(val_dict, err_val=-9999)
    mv.map_file(mukey_path, dst)


def mukey_rescue(gdb, raster_id=None, save=None,
                 exe_dir="~/github/ArcRasterRescue"):
    """Use the program Arc Raster Rescue from this fellow to extract the
    10m mukey raster from the ESRI Geodatabases (vectors are taken care of)
    
    https://github.com/r-barnes/ArcRasterRescue
    """

    # Expand user paths
    exe = os.path.join(os.path.expanduser(exe_dir), "arc_raster.exe")
    gdb = os.path.expanduser(gdb)
    save = os.path.expanduser(save)

    # Output file path
    if not gdb[-1] == "/":
        gdb = gdb + "/"

    # Try to detect rasters without raster_id, or rescue with
    if raster_id is None:
        call = [exe, gdb]
    else:
        if save is None:
            raise KeyError("No save path provided.")
        call = [exe, gdb, raster_id, save]

    try:
        x = sp.call(call)
    except Exception as e:
        print(e)
        print(RESCUE_INSTRUCTIONS)

    fix_mukey(save)



class Map_Soil:
    """Methods for mapping different soil properties at different depths.
    
    Sample Arguments
    ----------------
    mukey_path = "~/data/weto/soil/mukey_de.tif"
    gdb_path = "~/data/weto/soil/gNATSGO_DE.gdb"
    dst = "~/data/weto/soil/brockdepmin.tif"
    variable = "brockdepmin"
    """

    def __init__(self, gdb_path, mukey_path):
        """Initialize Map_Soil instance."""

        self.gdb_path = os.path.expanduser(gdb_path)
        self.mukey_path = os.path.expanduser(mukey_path)

    def __repr__(self):
        
        attrs = ["{}='{}'".format(k, v) for k, v in self.__dict__.items()]
        attrs_str = " ".join(attrs)
        msg = "<Map_Soil {}> ".format(attrs_str)

        return msg


    def set_table(self):
        """
        Create a raster attribute table with component and horizon values.
    
        Parameters
        ----------
        gdb_path : str
            Path to ESRI gSSRUGO Geodatabase
        mukey_path : str
            Path to a raster of the gSSURGO's 10 m map unit key raster.
        variable : str
            Name of the soil variable to be mapped.
        dst : str
            Path to destination file.
    
        Notes
        -----
        Variable descriptions:
            https://data.nal.usda.gov/system/files/SSURGO_Metadata_-_Table_Column_Descriptions.pdf#page=81
    
        Units:
            https://jneme910.github.io/CART/chapters/Soil_Propert_List_and_Definition
    

        """
    
        # Expand user path
        mukey_path = os.path.expanduser(mukey_path)
        gdb_path = os.path.expanduser(gdb_path)
        dst = os.path.expanduser(dst)
    
        # Get the Map Unit Aggregated Attribute Table
        mukey = xr.open_rasterio(mukey_path, chunks=(1, 5000, 5000))
        muaggatt = gpd.read_file(gdb_path, layer="muaggatt")
        chorizon = gpd.read_file(gdb_path, layer="chorizon")
        components = gpd.read_file(gdb_path, layer="component")
        components = pd.merge(chorizon, components, on="cokey")
        components = pd.merge(components, muaggatt, on="mukey")
    
        # Put the keys in front
        keys = [c for c in components.columns if "key" in c]
        others = [c for c in components.columns if "key" not in c]
        new_order = keys + others
        components = components[new_order]
    
        # Get the Horizon Table
        variable_df = components[["mukey", "chkey", "hzname", variable]]
        units = muaggatt[["mukey", "muname"]]
        variable_df = pd.merge(variable_df, units, on="mukey")
        variable_df = variable_df.dropna()
    
        # Now, whats the best way to map these values
        val_dict = dict(zip(variable_df["mukey"].astype(int),
                            variable_df[variable]))
        mv = Map_Values(val_dict, err_val=-9999)
        mv.map_file(mukey_path, dst)


class Site_Soil:
    '''
    Cortez Area Soil Survey:
    https://www.nrcs.usda.gov/Internet/FSE_MANUSCRIPTS/colorado/CO671/0/CO671%20Cortez.pdf

    SSURGO Metadata:
    http://www.nrcs.usda.gov/wps/PA_NRCSConsumption/download?cid=stelprdb1241115&ext=pdf
    https://data.nal.usda.gov/system/files/SSURGO_Metadata_-_Table_Column_Descriptions.pdf

    gSSURGO Source:
    https://nrcs.app.box.com/v/soils

    Minimum Data Needed:

        Upper and lower horizon depths (cm)
        Percentage sand, silt, and clay content
        1/3 bar bulk density
        Organic carbon
        pH in water
        Aluminum saturation
        Root abundance information

    Could also use:
        Soil series name
        Soil Classification
        Color
        Percent Slope
        Runoff Potential
        Fertility Factor
        Color
        Drainage
        Cation Exchange Capacity (cmol/kg)
        Percent Stones
        Percent Total Nitrogen
        Root Growth Factor (0 - 1)

    Okay, the current strategy is to download the gridded SSURGO (gSSURGO)
    geographic data base and work from the local file. To make this work, use
    ArcMap to pull in the 10m map from the geodatabase (One state at a time
    for now), join the components table with mukey and use that to join the
    chorizon table to that. This will give you all of the soil horizon needed
    to calculate the parameters needed for DSSAT, I think...hopefully.
    '''
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon
        self.soldf = pd.read_csv('data/co_soil.csv')
        self.solnc = xr.open_dataset('data/co_soil.nc')
        # self.densities = pd.read_csv('data/tables/bulk_densities.txt', sep='|')

    def nearest_albers(self, y, x):
        '''
        Take a coordinate from some where and return the index positions of the
        closest grid cell.
        '''
        data = self.solnc
        ys = data.y.values
        xs = data.x.values
        y_range = [ys[ys > y][0], ys[ys < y][-1]]
        y_diff = [abs(y-l) for l in y_range]
        y_idx = y_diff.index(min(y_diff))
        target_y = y_range[y_idx]
        x_range = [xs[xs > x][0], xs[xs < x][-1]]
        x_diff = [abs(x - l) for l in x_range]
        x_idx = x_diff.index(min(x_diff))
        target_x = x_range[x_idx]
        y = np.where(ys == target_y)[0][0]
        x = np.where(xs == target_x)[0][0]

        return y, x

    def get_mukey(self):
        """Get the MUKEY for a location."""

        # We have to project lat lon to albers equal area conic USGS
        wgs = Proj(init='epsg:4326')
        albers = Proj('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 ' +
                      '+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 ' +
                      '+towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
        x, y = transform(wgs, albers, self.lon, self.lat)
        y, x = self.nearest_albers(y, x)

        # Now find that spot in the soil netcdf
        data = self.solnc
        mukey = data.variables['Band1'][y, x].data

        return int(mukey)

    def get_horizon_depths(self):
        df = self.soldf
        mukey = self.get_mukey()
        hzdept_r = df['hzdept_r'][df['MUKEY'] == mukey].values
        hzdepb_r = df['hzdepb_r'][df['MUKEY'] == mukey].values

        return hzdept_r, hzdepb_r

    def get_bulk_density(self):
        df = self.soldf
        mukey = self.get_mukey()
        dbthirbar_r = df['dbthirdbar'][df['MUKEY'] == mukey].values
        return dbthirbar_r

    def get_percent_fines(self):
        df = self.soldf
        mukey = self.get_mukey()
        dbthirbar_r = self.get_bulk_density()
        sandtotal_r = df['sandtotal_'][df['MUKEY'] == mukey].values
        silttotal_r = df['silttotal_'][df['MUKEY'] == mukey].values
        claytotal_r = df['claytotal_'][df['MUKEY'] == mukey].values

        psand = sandtotal_r / dbthirbar_r  # Is this somehow already in %?
        psilt = silttotal_r / dbthirbar_r
        pclay = claytotal_r / dbthirbar_r

        return psand, psilt, pclay

    def get_map_unit(self):
        df = self.soldf
        mukey = self.get_mukey()
        muname = df['muname'][df['MUKEY'] == mukey].values

        return muname
