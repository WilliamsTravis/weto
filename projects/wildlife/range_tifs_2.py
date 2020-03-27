"""
Build overlapping wildlife range tiffs
"""

import os

from itertools import combinations

import dask.array as da
import numpy as np
import rasterio
from tqdm import tqdm
import xarray as xr

from numpy.random import randint
from dask.distributed import Client
from gdalmethods import Data_Path, to_raster


DP = Data_Path("/projects/rev/data/conus/wildlife")
DP2 = Data_Path("/scratch/twillia2/weto/wildlife")

RANGE_CATEGORIES = {
    1: {"eastern_red": "bat_ranges/eastern_red_bat_range.tif",
        "hoary": "bat_ranges/hoary_range.tif",
        "silver": "bat_ranges/silver_range.tif",
        # "big_brown": "bat_ranges/.tif",
        "tricolor": "bat_ranges/tricolored_bat_extent.tif",
        "little_brown": "bat_ranges/little_brown_bat_extent.tif",
        "brazilian": "bat_ranges/brazilian_free_tail_range.tif"},
    2: {"indiana": "bat_ranges/spext_indiana_bat.tif",
        "long_eared": "bat_ranges/northern_long_eared_bat_extent.tif",
        "gray": "bat_ranges/gray_myotis_range.tif"},
    3: {"condor": "bird_ranges/california_condor_range.tif",
        "crane": "bird_ranges/whooping_crane_range.tif",
        "owl": "bird_ranges/burrowing_owl_range.tif",
        "tortoise": "other_ranges/desert_tortoise_range.tif"},
    4: {"greater_sg": "bird_ranges/.greater_sg_rangetif",
        "greater_pc": "bird_ranges/greater_pc_range.tif",
        "lesser_pc": "bird_ranges/lesser_pc_range.tif",
        "gunnison_sg": "bird_ranges/gunnison_sg_range.tif"},
    5: {"golden": "bird_ranges/golden_eagle_range.tif",
        "bald": "bird_ranges/bald_eagle_range.tif"}
    }


CHUNKS = {"band": 1, "x": 5000, "y": 5000}




def combo_check(path_dict):
    """Check how many possible combinations of a set of ranges can be made.
    We can use this to check against the number of unique products we use
    as category values."""

    keys = path_dict.keys()
    
    combos = [[c for c in combinations(keys, i + 1)] for i in range(len(keys))]
    combos = [c for sc in combos for c in sc]
    combos = [sorted(c) for c in combos]
    combos = np.unique(combos)
    ncombos = len(combos)

    return ncombos


def code_range(path_dict):
    """
    path_dict = RANGE_CATEGORIES[5]
    """
    
    # Create a dictionary of random values 
    key_dict = {}
    for i, key in enumerate(path_dict.keys()):
        val = np.random.randint(1, 100)
        key_dict[key] = val
    try:
        assert len(path_dict) == len(key_dict)
    except AssertionError:
        msg = "\n\nNot all random keys are unique, try again."
        raise AssertionError(msg) 

    # We'll need to associate the numbers with categories
    number_dict = {k: i for i, k in key_dict.items()}
    vals = key_dict.values()
 
    # Find all unique products
    combos = [[c for c in combinations(vals, i + 1)] for i in range(len(vals))]
    combos = [c for sc in combos for c in sc]
    products = np.unique([np.prod(c) for c in combos])

    nproduct = len(products)
    check = combo_check(path_dict)
    try:
        assert nproduct == check
    except AssertionError:
        msg = "\n\nNot all products are unique, try something else."
        raise AssertionError(msg) 
    

    # combo_keys = {}
    # for combo in combos:
    #     key = "-".join([number_dict[c] for c in combo])
    #     value = seq_int(combo)
    #     combo_keys[key] = value

    # Assign each raster a unique value
    print("Stacking rasters and assigning unique values...")
    arrays = []
    for key, path in tqdm(path_dict.items(), position=0):
        value = key_dict[key]
        full_path = DP.join(path)
        ds = xr.open_rasterio(full_path)
        navalue = ds.attrs["nodatavals"][0]
        array = ds[0].data
        array[da.isnan(array)] = 0
        if not np.isnan(navalue):
            array[array == navalue] = 0
        array[array > 0] = value
        array[array == 0] = 1
        arrays.append(array)

    # Take the product, these will be in the code keys above
    stack = np.stack(arrays, axis=0)
    layer = np.prod(stack, axis=0)

    # These numbers are no good, we need to reset them
    uvals  = np.unique(layer)
    reset_dict = {val: i + 1  for i, val in enumerate(uvals)}
    for old, new in tqdm(reset_dict.items(), position=0):
        layer[layer == old] = new

    # These will now be small enough for this
    layer = layer.astype("uint8")

    # Save to temp and delete
    template = rasterio.open(full_path)
    temp_path = DP2.join("test.tif")
    with rasterio.Env():
        profile = template.profile
        profile["dtype"] = layer.dtype.name
        profile.update(
            dtype=rasterio.uint8,
            count=1,
            compress='lzw')
        with rasterio.open(temp_path, 'w', **profile) as dst:
            dst.write(layer, 1)


if __name__ == "__main__":
    path_dict = RANGE_CATEGORIES[5]
    code_range(path_dict)
