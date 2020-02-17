#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create masks for the ranking procedure. Separated here for memory purposes

Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""
import dask.array as da
import numpy as np
import os
from osgeo import gdal
from weto.functions import Data_Path, to_raster
import xarray as xr

# Move to functions
@da.as_gufunc(signature="(i)->(i)", output_dtypes=int, vectorize=True,
              allow_rechunk=True)
def gmask(array, navalues=[-9999., 0.0]):
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


# data path
dp = Data_Path("/scratch/twillia2/weto/data")

# Get each array, remove band dimension
blm = xr.open_rasterio(dp.join("rasters/albers/blm_codes_ac.tif"))[0]
tribes = xr.open_rasterio(dp.join("rasters/albers/tribal_codes_ac.tif"))[0]
state = xr.open_rasterio(dp.join("rasters/albers/state_land_codes_ac.tif"))[0]

# Make a mask of each higher priority layer
bmask = gmask(blm)
tmask = gmask(tribes)
smask = gmask(state)

# Make two composite masks, the max will bring 1s to the front
mask1 = da.stack([bmask, tmask, smask], axis=0).max(axis=0)
mask2 = da.stack([bmask, tmask], axis=0).max(axis=0)
mask3 = bmask

# Now we actually want the inverse of these (1 -> 0 -> 0 ; 0 -> -1 -> 1)
mask1 = (mask1 - 1) * -1
mask2 = (mask2 - 1) * -1
mask3 = (mask3 - 1) * -1

# Compute each of these and save to file
os.makedirs(dp.join("rasters/albers/masks/"), exist_ok=True)
reference = gdal.Open(dp.join("rasters/albers/blm_codes_ac.tif"))
proj = reference.GetProjection()
geom = reference.GetGeoTransform()

mask1 = mask1.compute().value.data
to_raster(mask1, dp.join("rasters/albers/masks/mask1.tif"), proj, geom)
del mask1

mask2 = mask2.compute().value.data
to_raster(mask2, dp.join("rasters/albers/masks/mask2.tif"), proj, geom)
del mask2

mask3 = mask3.compute()
to_raster(mask3, dp.join("rasters/albers/masks/mask3.tif"), proj, geom)
del mask3
