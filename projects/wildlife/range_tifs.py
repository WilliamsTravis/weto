"""
Build overlapping wildlife range tiffs
"""

import dask.array as da
from itertools import combinations
import numpy as np
import os
import rasterio
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


def seq_int(seq):
    """Concatenate a list of integers.
    seq = (5, 1, 4, 3)
    """
    if len(seq) > 0:
        if len(seq) > 1:
            frame = "".join(["{}" for i in range(len(seq))])
            return int(frame.format(*seq))
        else:
            return int("{}".format(seq[0])) 


def array_int(array):
    return np.apply_along_axis(seq_int, 0, array)


def code_range(path_dict):
    """
    path_dict = RANGE_CATEGORIES[5]
    """
    
    # Create a numeric dictionary for these keys
    key_dict = {key: i + 1  for i, key in enumerate(path_dict.keys())}
    number_dict = {k: i for i, k in key_dict.items()}
    vals = key_dict.values()
    combos = [[c for c in combinations(vals, i + 1)] for i in range(len(vals))]
    combos = [c for sc in combos for c in sc]
    combo_keys = {}
    for combo in combos:
        key = "-".join([number_dict[c] for c in combo])
        value = seq_int(combo)
        combo_keys[key] = value

    # Assign each raster a unique value
    arrays = []
    for key, path in path_dict.items():
        value = key_dict[key]
        full_path = DP.join(path)
        array = xr.open_rasterio(full_path, chunks=CHUNKS)[0].data
        array[da.isnan(array)] = 0
        array[array > 0] = value
        arrays.append(array)

    # Stack everything together - we might have to save this a temporary file
    stack = da.stack(arrays, axis=0)
    stack = stack.rechunk((stack.shape[0], 5000, 5000))
    stack = stack[:, 4000: 10000, 4000: 10000]

    # Try to map the function to each point
    client = Client()
    codes = da.apply_along_axis(seq_int, 0, stack, dtype="uint8")
    future = client.compute(codes)
    result = future.result()
    client.shutdown()
    client.close()

    # Save to temp and delete
    template = rasterio.open(full_path)
    temp_path = DP2.join("test.tif")
    with rasterio.Env():
        profile = template.profile
        profile.update(
            dtype=rasterio.uint8,
            count=1,
            compress='lzw')
        with rasterio.open(temp_path, 'w', **profile) as dst:
            dst.write(result)


if __name__ == "__main__":
    path_dict = RANGE_CATEGORIES[5]
    code_range(path_dict)
