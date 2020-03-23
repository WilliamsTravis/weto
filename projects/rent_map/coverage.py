#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate percent coverage of the cost code map.

Created on Thu Feb 20 22:26:30 2020

@author: twillia2
"""

import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr

from dask.distributed import Client
from gdalmethods import Data_Path


DP = Data_Path("~/data/weto/rent_map")


def codes():
    """Collect the cost codes from each category into a dictionary."""

    # Code to cost lookup table
    lookup = pd.read_csv(DP.join("tables/conus_cbe_lookup.csv"))

    # Split lookup table up into these categories
    blm = lookup[lookup["type"].str.contains("BLM Zone")]
    state = lookup[lookup["type"].str.contains("State Land")]
    private = lookup[(~lookup["type"].str.contains("BLM Zone")) &
                     (~lookup["type"].str.contains("Tribal Land")) &
                     (~lookup["type"].str.contains("State Land"))]
    tribal = lookup[lookup["type"].str.contains("Tribal Land")]

    # Assign each their own entry
    codes = {}
    codes["blm"] = blm["code"].values
    codes["state"] = state["code"].values
    codes["private"] = private["code"].values
    codes["tribal"] = tribal["code"].values

    return codes


def get_counts(codes):
    """Get cell counts for each category."""

    # Read in code and conus rasters
    chunks = {"band": 1, "x": 5000, "y": 5000}
    code_path = DP.join("rasters/albers/acre/cost_codes.tif")
    conus_path = DP.join("rasters/albers/acre/masks/conus.tif")
    codes = xr.open_rasterio(code_path, chunks=chunks)[0].data
    conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data
    coverage = codes[~np.isnan(codes)]
    conus = conus[~np.isnan(conus)]

    # Extract code from dictionary
    blm_codes = codes["blm"]
    tribal_codes = codes["tribal"]
    state_codes = codes["state"]
    private_codes = codes["private"]

    # Collect counts
    counts = {}
    with Client():

        # Arrays
        ntotal = da.count_nonzero(conus).compute()
        ncovered = da.count_nonzero(coverage).compute()
        excl_covered = coverage[coverage == 9999]
        blm_covered = coverage[da.isin(coverage, blm_codes)]
        tribal_covered = coverage[da.isin(coverage, tribal_codes)]
        state_covered = coverage[da.isin(coverage, state_codes)]
        private_covered = coverage[da.isin(coverage, private_codes)]

        # Counts
        counts["ntotal"] = ntotal
        counts["ncovered"] = ncovered
        counts["nexcl"] = da.count_nonzero(excl_covered).compute()
        counts["nexclnblm"] = da.count_nonzero(blm_covered).compute()
        counts["nexclntribal"] = da.count_nonzero(tribal_covered).compute()
        counts["nexclnstate"] = da.count_nonzero(state_covered).compute()
        counts["nexclnprivate"] = da.count_nonzero(private_covered).compute()

    return counts


def overall_coverage(counts):
    """Calculate percent coverage of each category overall."""
    
    # Overall percentages
    total_coverage = ncovered / ntotal
    blm_coverage = nblm / ntotal
    tribal_coverage = ntribal / ntotal
    private_coverage = nprivate / ntotal
    state_coverage = nstate / ntotal
    excl_coverage = nexcl / ntotal

    # Overall Data Frame
    overall = {"total": total_coverage,
               "blm": blm_coverage,
               "tribal": tribal_coverage,
               "state": state_coverage,
               "private": private_coverage,
               "excluded": excl_coverage}
    odf = pd.DataFrame(overall, index=[0]).T
    odf.columns = ["overall_percentage"]
    

def developable_coverage(counts):
    """Calculate percent coverage, but only within developable land."""

    # Developable Coverage
    covered = codes[(~np.isnan(codes)) & (codes != 9999)]
    developable = codes[(~np.isnan(codes)) & (codes != 9999)]
    developable[developable == 0] = 1
    with Client():
        ncovered = da.count_nonzero(covered).compute()
        ndevelopable = da.count_nonzero(developable).compute()

    # Developable Percentages
    total_coverage = ncovered / ndevelopable
    blm_coverage = nblm / ndevelopable
    tribal_coverage = ntribal / ndevelopable
    private_coverage = nprivate / ndevelopable
    state_coverage = nstate / ndevelopable

    # Developable Data Frame
    developable = {"total": total_coverage,
                   "blm": blm_coverage,
                   "tribal": tribal_coverage,
                   "state": state_coverage,
                   "private": private_coverage}
    ddf = pd.DataFrame(developable, index=[0]).T
    ddf.columns = ["developable_percentage"]

    return ddf


def coverage(counts):
    """Combine overall and developable coverage statistics."""

    # Get both data frames
    odf = overall_coverage(counts)
    ddf = developable_coverage(counts)

    # Composite Data Frame
    df = odf.join(ddf)
    df["area"] = df.index
    df = df[["area", "overall_percentage", "developable_percentage"]]
    df.to_csv(DP.join("tables/coverage.csv"), index=False)

    return df


if __name__ == "__main__":
    counts = get_counts()
    covdf = coverage(counts)
