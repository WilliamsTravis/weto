#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polygonizing the NLCD raster is too slow.

Created on Fri Jan 24 10:38:28 2020

@author: twillia2
"""
import subprocess as sp
import os
import rasterio
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
import sys


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


def polygonize(arg):
    """Use gdal to cut a raster into a smaller piece and polygonize
    
    Too noisy.
    
    extent = extents[0]
    """
    os.makedirs("chunks", exist_ok=True)

    # Separate arguments
    extent = arg[0]
    rfile = arg[1]
    chunk = arg[2]

    # Get everything in order
    extent = [str(e) for e in extent]
    chunk = "{:02d}".format(chunk)
    folder = os.getcwd()
    routfile = os.path.basename(rfile).split(".")[0] + "_" + chunk + ".tif"
    poutfile = os.path.basename(rfile).split(".")[0] + "_" + chunk + ".shp"
    rout = os.path.join(folder, "chunks", routfile)  
    pout = os.path.join(folder, "chunks", poutfile)

    # Let's not overwrite
    if not os.path.exists(rout):
        sp.call(["gdalwarp", "-q", "-te", extent[0], extent[1], extent[2],
                 extent[3], rfile, rout], stdout=sp.PIPE, stderr=sp.PIPE)

    if not os.path.exists(pout):
        sp.call(["gdal_polygonize.py", "-8", "-q", rout, pout],
                stdout=sp.PIPE, stderr=sp.PIPE)

    return pout


def main(extents, rfile):
    ncpu = mp.cpu_count()
    chunknumbers = [i for i in range(len(extents))]
    rfiles = np.repeat(rfile, len(extents))
    args = list(zip(extents, rfiles, chunknumbers))
    with mp.Pool(ncpu) as pool:
        pfiles = []
        for pfile in tqdm(pool.imap(polygonize, args),
                      total=len(extents[:10]), position=0, file=sys.stdout):
            pfiles.append(pfile)

    return pfiles




# To do:

# sp.call(["ogr2ogr", "-f", "ESRI Shapefile", "-update", "-append",
#          "nlcd_2018_ag.shp", "chunks/nlcd_2016_ag_*shp", "-nln", "DN"])

# sp.call(['ogr2ogr', 'nlcd_2018_ag_dissolved.shp', 'nlcd_2018_ag.shp',
#          '-dialect', 'sqlite', '-sql', '"SELECT ST_Union(geometry), 
#           DN FROM input GROUP BY DN"'])


if __name__ == '__main__':
    os.chdir(".")
    RFILE = "/scratch/twillia2/weto_costs/data/rasters/nlcd_2016_ag.tif"
    print("Running nlcd_polygonize on " + RFILE)
    sp.call(["gdalinfo", "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
    extents = split_extent(RFILE, n=225)
    main(extents, RFILE)
