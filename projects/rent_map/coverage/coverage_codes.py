# -*- coding: utf-8 -*-
"""
Calculate coverage statistics for codes.

Created on Thu Feb 20 22:26:30 2020

@author: twillia2
"""

import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path
from tqdm import tqdm


# DP = Data_Path("~/data/weto/rent_map")
DP = Data_Path("/scratch/twillia2/weto/data")


def fixit(x):
    x = x.replace("*", "").replace("\n", "").replace("$", "").replace(",", "")
    if x == " NULL " or  x == "  NULL  ":
        x = 0
    return x


def get_codes(cost_coverage=False):
    """Collect the cost codes from each category into a dictionary."""

    # Code to cost lookup table
    lookup = pd.read_csv(DP.join("tables/conus_cbe_lookup.csv"))
    lookup.columns = ['code', 'type', 'dollar_ac']
    lookup["dollar_ac"] = lookup["dollar_ac"].apply(fixit)
    lookup["dollar_ac"] = lookup["dollar_ac"].astype(float)

    if cost_coverage:
        lookup = lookup[lookup["dollar_ac"] > 0.0]

    # Split lookup table up into these categories
    blm = lookup[lookup["type"].str.contains("BLM Zone")]
    state = lookup[lookup["type"].str.contains("State Land")]
    private = lookup[(~lookup["type"].str.contains("BLM Zone")) &
                      (~lookup["type"].str.contains("Tribal Land")) &
                      (~lookup["type"].str.contains("State Land"))]
    tribal = lookup[lookup["type"].str.contains("Tribal Land")]

    # Assign each their own entry
    code_dict = {}
    code_dict["blm"] = blm["code"].values
    code_dict["state"] = state["code"].values
    code_dict["private"] = private["code"].values
    code_dict["tribal"] = tribal["code"].values

    return code_dict


def get_counts(cost_coverage=False):
    """Get cell counts for each category."""

    code_dict = get_codes(cost_coverage)

    # Read in code and conus rasters
    chunks = {"band": 1, "x": 5000, "y": 5000}
    code_path = DP.join("rasters/albers/acre/cost_codes.tif")
    cost_path = DP.join("rasters/albers/acre/rent_map.tif")
    conus_path = DP.join("rasters/albers/acre/masks/conus.tif")
    codes = xr.open_rasterio(code_path, chunks=chunks)[0].data
    costs = xr.open_rasterio(cost_path, chunks=chunks)[0].data
    conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data

    # Dask array's `count_nonzero` counts na values
    codes[da.isnan(codes)] = 0
    conus[da.isnan(conus)] = 0

    # If calculating costs
    if cost_coverage:
        coverage = codes[(costs > 0) | (codes == 9999)]  # No exclusion in cost
    else:
        coverage = codes.copy()

    # Extract code from dictionary
    blm_codes = code_dict["blm"]
    tribal_codes = code_dict["tribal"]
    state_codes = code_dict["state"]
    private_codes = code_dict["private"]

    # Arrays
    developable = conus[codes != 9999]
    dev_covered = coverage[coverage != 9999]
    excl = coverage[coverage == 9999]
    blm = coverage[da.isin(coverage, blm_codes)]
    tribal = coverage[da.isin(coverage, tribal_codes)]
    state = coverage[da.isin(coverage, state_codes)]
    private = coverage[da.isin(coverage, private_codes)]
    arrays = {"excl": excl, "blm": blm, "tribal": tribal, "state": state,
              "private": private, "covered": coverage, "total": conus, 
              "developable": developable, "dev_covered": dev_covered}

    # Collect counts
    counts = {}
    with Client():
        for key, item in tqdm(arrays.items(), position=0):
            counts["n" + key] = da.count_nonzero(item).compute()

    return counts


def overall_coverage(counts):
    """Calculate percent coverage of each category overall."""

    # Overall percentages
    total_coverage = counts["ncovered"] / counts["ntotal"]
    blm_coverage = counts["nblm"] /counts["ntotal"]
    tribal_coverage = counts["ntribal"] / counts["ntotal"]
    private_coverage = counts["nprivate"] / counts["ntotal"]
    state_coverage = counts["nstate"] / counts["ntotal"]
    excl_coverage = counts["nexcl"] / counts["ntotal"]

    # Overall Data Frame
    overall = {"total": total_coverage,
                "blm": blm_coverage,
                "tribal": tribal_coverage,
                "state": state_coverage,
                "private": private_coverage,
                "excluded": excl_coverage}
    odf = pd.DataFrame(overall, index=[0]).T
    odf.columns = ["overall_percentage"]

    return odf


def developable_coverage(counts):
    """Calculate percent coverage, but only within developable land."""

    # Developable Percentages
    total_coverage = counts["ndev_covered"] / counts["ndevelopable"]
    blm_coverage = counts["nblm"] / counts["ndevelopable"]
    tribal_coverage = counts["ntribal"] / counts["ndevelopable"]
    private_coverage = counts["nprivate"] / counts["ndevelopable"]
    state_coverage = counts["nstate"] / counts["ndevelopable"]

    # Developable Data Frame
    developable = {"total": total_coverage,
                    "blm": blm_coverage,
                    "tribal": tribal_coverage,
                    "state": state_coverage,
                    "private": private_coverage}
    ddf = pd.DataFrame(developable, index=[0]).T
    ddf.columns = ["developable_percentage"]

    return ddf


def get_coverage(counts, file_path="coverage.csv"):
    """Combine overall and developable coverage statistics."""

    # Get both data frames
    odf = overall_coverage(counts)
    ddf = developable_coverage(counts)

    # Composite Data Frame
    df = odf.join(ddf)
    df["area"] = df.index
    df = df[["area", "overall_percentage", "developable_percentage"]]
    df.to_csv(DP.join("tables", file_path), index=False)

    return df


if __name__ == "__main__":
    counts = get_counts(cost_coverage=False)
    covdf = get_coverage(counts, file_path="coverage_codes.csv")
