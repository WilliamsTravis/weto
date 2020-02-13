#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:37:42 2020

@author: travis
"""

import dask.array as da
import rasterio
from dask import is_dask_collection
from dask.base import tokenize
from rasterio.windows import Window


# READ
def read_raster(path, band=None, block_size=1):
    """Read all or some bands from raster
    Arguments:
        path {string} -- path to raster file
    Keyword Arguments:
        band {int, iterable(int)} -- band number or iterable of bands.
            When passing None, it reads all bands (default: {None})
        block_size {int} -- block size multiplier (default: {1})
    Returns:
        dask.array.Array -- a Dask array
    """

    if isinstance(band, int):
        return read_raster_band(path, band=band, block_size=block_size)
    else:
        if band is None:
            bands = range(1, get_band_count(path) + 1)
        else:
            bands = list(band)
        return da.stack([
            read_raster_band(path, band=band, block_size=block_size)
            for band in bands
        ])


def read_raster_band(path, band=1, block_size=1):
    """Read a raster band and return a Dask array
    Arguments:
        path {string} -- path to the raster file
    Keyword Arguments:
        band {int} -- number of band to read (default: {1})
        block_size {int} -- block size multiplier (default: {1})
    """

    def read_window(raster_path, window, band):
        with rasterio.open(raster_path) as src:
            return src.read(band, window=window)

    def resize_window(window, block_size):
        return Window(
            col_off=window.col_off * block_size,
            row_off=window.row_off * block_size,
            width=window.width * block_size,
            height=window.height * block_size)

    def block_windows(dataset, band, block_size):
        return [(pos, resize_window(win, block_size))
                for pos, win in dataset.block_windows(band)]

    with rasterio.open(path) as src:
        h, w = src.block_shapes[band - 1]
        chunks = (h * block_size, w * block_size)
        name = 'raster-{}'.format(tokenize(path, band, chunks))
        dtype = src.dtypes[band - 1]
        shape = src.shape
        blocks = block_windows(src, band, block_size)

    dsk = {(name, i, j): (read_window, path, window, band)
           for (i, j), window in blocks}

    return da.Array(dsk, name, chunks, dtype, shape)


def get_band_count(raster_path):
    """Read raster band count"""
    with rasterio.open(raster_path) as src:
        return src.count


# WRITE
def write_raster(path, array, **kwargs):
    """Write a dask array to a raster file
    If array is 2d, write array on band 1.
    If array is 3d, write data on each band
    Arguments:
        path {string} -- path of raster to write
        array {dask.array.Array} -- band array
        kwargs {dict} -- keyword arguments to delegate to rasterio.open
    Examples:
        # Write a single band raster
        >> red_band = read_raster_band("test.tif", band=1)
        >> write_raster("new.tif", red_band)
        # Write a multiband raster
        >> img = read_raster("test.tif")
        >> new_img = process(img)
        >> write_raster("new.tif", new_img)
    """
    if len(array.shape) != 2 and len(array.shape) != 3:
        raise TypeError('invalid shape (must be either 2d or 3d)')

    if is_dask_collection(array):
        with RasterioDataset(path, 'w', **kwargs) as dst:
            da.store(array, dst, lock=True)
    else:
        with rasterio.open(path, 'w', **kwargs) as dst:
            if len(array.shape) == 2:
                dst.write(array, 1)
            else:
                dst.write(array)


class RasterioDataset:
    """Rasterio wrapper to allow dask.array.store to do window saving.
    Example:
        >> rows = cols = 21696
        >> a = da.ones((4, rows, cols), dtype=np.float64, chunks=(1, 4096, 4096) )
        >> a = a * np.array([255., 255., 255., 255.])[:, None, None]
        >> a = a.astype(np.uint8)
        >> with RasterioDataset('test.tif', 'w', driver='GTiff', width=cols, height=rows, count=4, dtype=np.uint8) as r_file:
        ..    da.store(a, r_file, lock=True)
    """

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.dataset = None

    def __setitem__(self, key, item):
        """Put the data chunk in the image"""
        if len(key) == 3:
            index_range, y, x = key
            indexes = list(
                range(index_range.start + 1, index_range.stop + 1,
                      index_range.step or 1))
        else:
            indexes = 1
            y, x = key

        chy_off = y.start
        chy = y.stop - y.start
        chx_off = x.start
        chx = x.stop - x.start

        self.dataset.write(
            item, window=Window(chx_off, chy_off, chx, chy), indexes=indexes)

    def __enter__(self):
        """Enter method"""
        self.dataset = rasterio.open(*self.args, **self.kwargs)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit method"""
        self.dataset.close()