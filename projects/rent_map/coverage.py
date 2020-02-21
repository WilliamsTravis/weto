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

dp = Data_Path("/scratch/twillia2/weto/data")
code_path = dp.join("rasters/albers/acre/cost_codes.tif")
conus_path = dp.join("rasters/albers/acre/masks/conus.tif")
chunks = {"band": 1, "x": 5000, "y": 5000}
codes = xr.open_rasterio(code_path, chunks=chunks)[0].data
conus = xr.open_rasterio(conus_path, chunks=chunks)[0].data

# For code catgories
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
blm = lookup[lookup["type"].str.contains("BLM Zone")]
state = lookup[lookup["type"].str.contains("State Land")]
private = lookup[(~lookup["type"].str.contains("BLM Zone")) &
                 (~lookup["type"].str.contains("Tribal Land")) &
                 (~lookup["type"].str.contains("State Land"))]
tribal = lookup[lookup["type"].str.contains("Tribal Land")]
blm_codes = blm["code"].values
state_codes = state["code"].values
private_codes = private["code"].values
tribal_codes = tribal["code"].values


# Overall counts
coverage = codes[~np.isnan(codes)]
conus = conus[~np.isnan(conus)]
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
    nexcl = da.count_nonzero(excl_covered).compute()
    nblm = da.count_nonzero(blm_covered).compute()
    ntribal = da.count_nonzero(tribal_covered).compute()
    nstate = da.count_nonzero(state_covered).compute()
    nprivate = da.count_nonzero(private_covered).compute()

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
           "excluded": excl_coverage
           }
odf = pd.DataFrame(overall, index=[0]).T
odf.columns = ["overall_percentage"]

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
               "private": private_coverage
               }
ddf = pd.DataFrame(developable, index=[0]).T
ddf.columns = ["developable_percentage"]

# Composite Data Frame
df = odf.join(ddf)
df.to_csv(dp.join("tables/coverage.csv"), index=False)
