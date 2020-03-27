#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate coverage by census divisions.

Created on Tue Mar 24 10:23:16 2020

@author: twillia2
"""

import os

from tqdm import tqdm

import dask.array as da
import geopandas as gpd
import pandas as pd
import rasterio
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path, rasterize


DP = Data_Path("/scratch/twillia2/weto/data")
TEMPLATE = DP.join("rasters", "albers", "acre", "rent_map.tif")
DIVISIONS_PATH = DP.join("rasters", "albers", "acre", "census_divisions.tif")


# Get a shapefile of the census blocks
def make_divisions():
    """Download reproject and rasterize the US Census divisions."""

    # Paths
    shp_path = DP.join("shapefiles", "USA", "cb_2018_us_division_5m.shp")
    tproj = rasterio.open(TEMPLATE).crs.to_proj4()

    if not os.path.exists(shp_path):

        # Read
        print("Retrieving Census Divisions...")
        divisions = gpd.read_file("https://www2.census.gov/geo/tiger/"
                                  "GENZ2018/shp/cb_2018_us_division_5m.zip")

        # Reproject
        divisions = divisions.to_crs(tproj)

        # Save
        divisions.to_file(shp_path)

    if not os.path.exists(DIVISIONS_PATH):

        # Rasterize
        print("Rasterizing Census Divisions...")
        rasterize(shp_path, DIVISIONS_PATH, attribute="DIVISIONCE",
                  template_path=TEMPLATE)

    # Associate division codes with descriptions
    divisions = gpd.read_file(shp_path)
    divisions["DIVISIONCE"] = divisions["DIVISIONCE"].astype(int)
    division_dict = dict(zip(divisions["DIVISIONCE"], divisions["NAME"]))

    return division_dict


def coverage_total(division_dict, cost_coverage=False):

    conus_path = DP.join("rasters", "albers", "acre", "masks", "conus.tif")
    code_path = DP.join("rasters", "albers", "acre", "cost_codes.tif")
    cost_path = DP.join("rasters/albers/acre/rent_map.tif")
    chunks = {"band": 1, "x": 5000, "y": 5000}

    # Read in the tifs
    codes = xr.open_rasterio(code_path, chunks=chunks)[0].data
    conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data
    costs = xr.open_rasterio(cost_path, chunks=chunks)[0].data
    divisions = xr.open_rasterio(DIVISIONS_PATH, chunks=chunks)[0].data

    # Set nans to zero (count_nonzero counts nans)
    codes[da.isnan(codes)] = 0

    coverages = {}
    with Client():
        for key, item in tqdm(division_dict.items(), position=0):
            
            div = conus[divisions == key]
            total = da.count_nonzero(div)

            # If calculating costs
            if cost_coverage:
                coverage = codes[((costs > 0) | (codes == 9999)) &
                                 (divisions == key)]
                coded = da.count_nonzero(coverage)
            else:
                coverage = codes[divisions == key]
                coded = da.count_nonzero(coverage)
            ratio = coded / total
            coverages[item] = ratio.compute()

    df = pd.DataFrame(coverages, index=[0]).T
    df.columns = ["total_coverage"]

    return df


def coverage_developable(division_dict, cost_coverage=False):

    conus_path = DP.join("rasters", "albers", "acre", "masks", "conus.tif")
    code_path = DP.join("rasters", "albers", "acre", "cost_codes.tif")
    cost_path = DP.join("rasters/albers/acre/rent_map.tif")
    chunks = {"band": 1, "x": 5000, "y": 5000}

    # Read in the tifs
    codes = xr.open_rasterio(code_path, chunks=chunks)[0].data
    conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data
    costs = xr.open_rasterio(cost_path, chunks=chunks)[0].data
    divisions = xr.open_rasterio(DIVISIONS_PATH, chunks=chunks)[0].data

    # Set nans to zero (count_nonzero counts nans)
    codes[da.isnan(codes)] = 0

    coverages = {}
    with Client():
        for key, item in tqdm(division_dict.items(), position=0):
            # break
            developable = da.count_nonzero(
                conus[(divisions == key) &
                      (codes != 9999)]
                )

            # If calculating costs
            if cost_coverage:
                coverage = codes[(costs > 0) &
                                 (codes != 9999) &
                                 (divisions == key)]
                coded = da.count_nonzero(coverage)
            else:
                coverage = codes[(divisions == key) &
                                 (codes != 9999)]
                coded = da.count_nonzero(coverage)
            ratio = coded / developable
            coverages[item] = ratio.compute()

    df = pd.DataFrame(coverages, index=[0]).T
    df.columns = ["developable_coverage"]

    return df

def merge_dfs(df1, df2):

    df = df1.join(df2)
    df["census_division"] = df.index
    cols = df.columns
    new_cols = [cols[i] for i in (2, 0, 1)]
    df = df[new_cols]
    
    return df


if __name__ == "__main__":

    division_dict = make_divisions()

    # Codes
    tcode_df = coverage_total(division_dict, cost_coverage=False)
    dcode_df = coverage_developable(division_dict, cost_coverage=False)
    codedf = merge_dfs(tcode_df, dcode_df)
    save_path = DP.join("tables", "census_code_coverage.csv")
    codedf.to_csv(save_path, index=False)

    # Costs
    tcost_df = coverage_total(division_dict, cost_coverage=True)
    dcost_df = coverage_developable(division_dict, cost_coverage=True)
    costdf = merge_dfs(tcost_df, dcost_df)
    save_path = DP.join("tables", "census_cost_coverage.csv")
    costdf.to_csv(save_path, index=False)
