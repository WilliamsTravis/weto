#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Manipulating the rasterstats author's own example mulitprocessing script. 

https://github.com/perrygeo/python-rasterstats/blob/master/examples/multiproc.py

Created on Thu Feb 27 12:50:56 2020

@author: twillia2
"""

import fiona
import itertools
import multiprocessing as mp
import numpy as np
import rasterio
import sys
import p_tqdm

from gdalmethods import Data_Path
from rasterstats import zonal_stats
from tqdm import tqdm


# Paths
DPA = Data_Path("/scratch/twillia2/weto/populations/albers")
TIF = DPA.join("landscan_night_2017.tif")
SHP = DPA.join("sc_circle_single_buffer.gpkg")
NA = rasterio.open(TIF).nodata


def get_chunks(data, n):
    """Yield successive n-sized chunks from a slice-able iterable."""
    chunks = []
    for i in range(0, len(data), n):
        chunks.append(data[i: i + n])

    return chunks


def zonal_stats_partial(feats):
    """Wrapper for zonal stats, takes a list of features"""

    stats = zonal_stats(feats, TIF, nodata=NA, all_touch=True,
                        stats="sum max min")
    out = [{**feats[i], **stats[i]} for i, _ in enumerate(feats)]

    return out


if __name__ == "__main__":

    # Create a process pool using all cores minus 1  # <----------------------- How to manage this sort of thing with SLURM?
    ncpu = mp.cpu_count() - 1

    # Open list of features
    with fiona.open(SHP) as src:
        features = list(src)

    # Parallel
    chunks = get_chunks(features, ncpu)[:3]
    # stats_lists = []
    # with mp.Pool(ncpu) as pool:
    #     for out in tqdm(pool.imap(zonal_stats_partial, chunks), position=0,
    #                     total=len(chunks), file=sys.stdout):
    #         stats_lists.append(out)   

    # pool = mp.Pool(ncpu)
    # stats_list = pool.map(zonal_stats_partial, chunks)

    # Trying a new package
    stats_lists = p_tqdm.p_umap(zonal_stats_partial, chunks)




    # # Serial
    # for f in tqdm(features):
    #     stats_lists.append(zonal_stats(f, TIF))

    # flatten to a single list
    stats = list(itertools.chain(*stats_lists))

    try:
        assert len(stats) == len(features)
        # save to file?
    except AssertionError:
        print("Error: not all features completed.")
