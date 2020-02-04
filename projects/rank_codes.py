#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a single field of lookup codes based on a ranking system.

1) read all four layers in as arrays
2) write a function that selects the highest ranked of those present at any
    single cell.
3) write to raster! That's the plan at least'

blm_codes_ac.tif   blm_codes.tif             nlcd_codes_ac.tif   state_land_codes_ac.tif   tribal_codes_ac.tif


Created on Mon Feb  3 16:02:50 2020

@author: twillia2
"""
import numpy as np
import rasterio
from weto.functions import Data_Path, to_raster, read_raster

# data path
dp = Data_Path("/scratch/twillia2/weto/data")

# best way to do this?
#1) Apply function to each point (computationally intensive)

#2) Apply cumulative masks (trickier, but less computation)
# Get each array
blm = rasterio.open(dp.join("rasters/albers/blm_codes_ac.tif")).read()
tribes = rasterio.open(dp.join("rasters/albers/tribal_codes_ac.tif")).read()
state = rasterio.open(dp.join("rasters/albers/state_land_codes_ac.tif")).read()
nlcd, geom, proj = read_raster(dp.join("rasters/albers/nlcd_codes_ac.tif"))

# Create a mask out of each
bmask = blm.copy()
bmask[bmask == -9999] = np.nan
bmask = bmask * 0 + 1

tmask = tribes.copy()
tmask[tmask == -9999] = np.nan
tmask = tmask * 0 + 1

smask = state.copy()
smask[smask == -9999] = np.nan
smask = smask * 0 + 1

nmask = nlcd.copy()
nmask[nmask == -9999] = np.nan
nmask = nmask * 0 + 1

# Then return nans to 0

# say we had three masks
# mask 1 = bmask + tmask + smask
# mask 2 = bmask + tmask 
# mask 3 = bmask

# Then we can iterate like this
# nlcd = nlcd * mask 1
# state = state * mask 2
# tribes = tribes * mask 3
# blm = blm

# And then fill in the gaps



