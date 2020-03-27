# -*- coding: utf-8 -*-
"""
Create the cost_codes.tif data set.

Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""

import os

import dask.array as da
import numpy as np
import xarray as xr

from osgeo import gdal

from dask.distributed import Client
from gdalmethods import Data_Path, to_raster
from tqdm import tqdm

# Data Paths
DP = Data_Path("/scratch/twillia2/weto/data")
DPM = Data_Path(DP.join("rasters/albers/acre/masks"))

# Set the chunk sizes
CHUNKS = {'band': 1, 'x': 5000, 'y': 5000}

# Reference geometry
REFERENCE = gdal.Open(DP.join("rasters/albers/acre/blm_codes.tif"))
PROJ = REFERENCE.GetProjection()
GEOM = REFERENCE.GetGeoTransform()


@da.as_gufunc(signature="(i)->(i)", output_dtypes=int, vectorize=True,
              allow_rechunk=True)
def gmask(array, navalues=(-9999., 0.0)):
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


# Make a conus mask with 1s for in and nans for out
def conus_mask():
    if not os.path.exists(DPM.join("conus.tif")):
        print("Creating CONUS mask...")
        template = DP.join("rasters/albers/acre/nlcd.tif")
        conus = xr.open_rasterio(template, chunks=CHUNKS)[0].data
        conus = conus.astype("float32")
        conus[conus == 0.] = np.nan
        conus = (conus * 0) + 1
        with Client():
            conus = conus.compute()
        to_raster(conus, DPM.join("conus.tif"), PROJ, GEOM, dtype="float32")


# Code paths
def code_paths():
    paths = {}
    paths["nlcd_path"] = DP.join("rasters/albers/acre/nlcd_codes.tif")
    paths["blm_path"] = DP.join("rasters/albers/acre/blm_codes.tif")
    paths["tribal_path"] = DP.join("rasters/albers/acre/tribal_codes.tif")
    paths["state_path"] = DP.join("rasters/albers/acre/state_codes.tif")
    paths["excl_path"] = DP.join("rasters/albers/acre/rent_exclusions.tif")

    return paths


def build_masks():

    # Pull paths out
    paths = code_paths()
    nlcd_path = paths["nlcd_path"]
    blm_path = paths["blm_path"]
    tribal_path = paths["tribal_path"]
    state_path = paths["state_path"]
    excl_path = paths["excl_path"]

    # Get each array, remove band dimension
    nlcd = xr.open_rasterio(nlcd_path, chunks=CHUNKS)[0].data
    blm = xr.open_rasterio(blm_path, chunks=CHUNKS)[0].data
    tribes = xr.open_rasterio(tribal_path, chunks=CHUNKS)[0].data
    state = xr.open_rasterio(state_path, chunks=CHUNKS)[0].data
    excl = xr.open_rasterio(excl_path, chunks=CHUNKS)[0].data

    # I need to fix this...nlcd's shape is off by one cell
    if nlcd.shape != blm.shape:
        # One too many on top
        nlcd = nlcd[1:, :]

        # Now we need two on bottom
        x = da.from_array(np.zeros((1, nlcd.shape[1])))
        nlcd = np.vstack((nlcd, x, x))

        # and one extra on the side
        y = da.from_array(np.zeros((nlcd.shape[0], 1)))
        nlcd = np.hstack((nlcd, y))
        with Client():
            nlcd = nlcd.compute()
        to_raster(nlcd, nlcd_path, PROJ, GEOM)

    # Make a mask of each higher priority layer
    emask = gmask(excl)
    bmask = gmask(blm)
    tmask = gmask(tribes)
    smask = gmask(state)

    # Make four composite masks, the max will bring 1s to the front
    mask1 = da.stack([emask, bmask, tmask, smask], axis=0).max(axis=0)
    mask2 = da.stack([emask, bmask, tmask], axis=0).max(axis=0)
    mask3 = da.stack([emask, bmask], axis=0).max(axis=0)
    mask4 = emask

    # Now we actually want the inverse of these (1 -> 0 -> 0 ; 0 -> -1 -> 1)
    mask1 = (mask1 - 1) * -1
    mask2 = (mask2 - 1) * -1
    mask3 = (mask3 - 1) * -1
    mask4 = (mask4 - 1) * -1

    # We'll need to setup interim files, I can't handle the memory useage
    layers = {"nlcd": nlcd,
              "state": state,
              "tribes": tribes,
              "blm": blm}
    masked_paths = {"nlcd": DPM.join("nlcd.tif"),
                    "state": DPM.join("state.tif"),
                    "tribes": DPM.join("tribes.tif"),
                    "blm": DPM.join("blm.tif")}
    masks = {"nlcd": mask1,
             "state": mask2,
             "tribes": mask3,
             "blm": mask4}

    # Loop through, the lowest priority layers are masked by the largest mask <---- This is obviously not the best way to do this
    print("Masking individual layers...")
    for key, layer in tqdm(layers.items(), position=0):
        layer_path = masked_paths[key]
        mask = masks[key]
        mlayer = layer * mask
        with Client():
            masked_layer = mlayer.compute()
        to_raster(masked_layer, layer_path, PROJ, GEOM)
        del masked_layer

    return masked_paths


def composite(masked_paths):
    
    # Too much at once, read the precomputd masked layers from file
    excl_path = DP.join("rasters/albers/acre/rent_exclusions.tif")
    excl = xr.open_rasterio(excl_path, chunks=CHUNKS)[0].data
    mblm = xr.open_rasterio(masked_paths["blm"], chunks=CHUNKS)[0].data
    mtribes = xr.open_rasterio(masked_paths["tribes"], chunks=CHUNKS)[0].data
    mstate = xr.open_rasterio(masked_paths["state"], chunks=CHUNKS)[0].data
    mnlcd = xr.open_rasterio(masked_paths["nlcd"], chunks=CHUNKS)[0].data
    mconus = xr.open_rasterio(DPM.join("conus.tif"), chunks=CHUNKS)[0].data

    print("Merging layers ...")
    with Client():
        layers = [excl, mblm, mtribes, mstate, mnlcd]
        composite_layer = da.stack(layers, axis=0).max(axis=0)
        composite_layer = composite_layer * mconus
        final = composite_layer.compute()

    # Save to file
    save = DP.join("rasters/albers/acre/cost_codes.tif")
    print("Saving file to " + save + " ...")
    to_raster(final, save, PROJ, GEOM)


def main():
    conus_mask()
    masked_paths = build_masks()
    composite(masked_paths)


if __name__ == "__main__":
    main()
