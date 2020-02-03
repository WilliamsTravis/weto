
# SBATCH --time=01:00:00
# SBATCH --acount=rev
# SBATCH --job-name=nlcd_map 
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=1

import numpy as np
import pandas as pd
import xarray as xr
from weto.functions import Data_Path

# Get the raster and the raster lookup values
dp = Data_Path("/scratch/twillia2/weto/data")
reference = pd.read_csv(dp.join("tables/nlcd_rast_lookup.csv"))
agcounty = xr.open_rasterio(dp.join("rasters/agcounty_product.tif"))
agcounty = agcounty[0, 40000:45000, 40000:45000]  # testing
map_vals = dict(zip(reference["rast_val"], reference["index"]))
map_vals[0] = 0

# Try this guy out
try:
    out = np.vectorize(map_vals.get)(agcounty)
    print("success.")
except Exception as e:
    print(e)
