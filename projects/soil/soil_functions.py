#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Downloading and fixing the Argonne wind mapper data.

Created on Mon Mar 16 13:40:17 2020

@author: travis
"""

import multiprocessing as mp
import os
import sys
import urllib.request
import zipfile as zp

from glob import glob
from itertools import product
from tqdm import tqdm

import fiona
import geopandas as gpd

from gdalmethods import Data_Path


DP = Data_Path("/home/travis/data/weto/blm_siting")
SRCDIR = DP.join("WWMP_Data_ViewerV2")
TRGTDIR = DP.join("WWMP_Data_ViewerV2/gpkgs")


def asc_converter(gdb, target_folder, overwrite=True):
    """Convert all datasets in an ESRI geogdatabase file to a geopackage and
    assign proper CRS."""

    # Make target folder
    trgt = Data_Path(target_folder)

    # Get the subdatasets
    print("Converting layers in " + gdb)
    layerlist = fiona.listlayers(gdb)
    for layer in tqdm(layerlist, position=0, file=sys.stdout):
        path = trgt.join(layer + ".gpkg")
        if not os.path.exists(path) and not overwrite:
            print(path + " exists. Use 'overwrite=True'")
            return

        data = gpd.read_file(gdb, layer=layer)
        data = data.to_crs("epsg:4326")
        data.to_file(path, driver="GPKG")


def get_wwmp():
    """Download and unzip the Argonne National Lab's West Wide Wind Mapping
    Project Data."""

    url = "http://wwmp.anl.gov/downloads/gis/WWMP_Data_ViewerV2.zip"
    dst = DP.join(os.path.basename(url))

    if not os.path.exists(dst):
        urllib.request.urlretrieve(url, dst)
        with zp.ZipFile(dst, "r") as zipf:
            zipf.extractall(os.path.dirname(dst))


def convert_all(src_folder, trgt_folder):
    """Convert all ESRI geodatabases in a folder to geopackages in another."""

    # Get geodatabase paths
    gdbs = glob(os.path.join(src_folder, "*gdb"))

    # Convert
    with mp.Pool(mp.cpu_count() - 1) as pool:
        pool.starmap(asc_converter, product(gdbs, trgt_folder))


if __name__ == "__main__":
    convert_all(SRCDIR, TRGTDIR)
