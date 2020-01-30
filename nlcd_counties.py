"""
Create a US cost surface layer based on NLCD land use type and county.

1) Get 2016 NLCD raster.
    - s3-us-west-2.amazonaws.com/mrlc/NLCD_2016_Land_Cover_L48_20190424.zip
    - 30 meter resolution.
    - This has a unique Albers Concical Equal Area projection and under
      referenced. Try projecting to WGS 84 with nearest neighbor to maintain
      NLCD category values and match other data sets. 
2) Create a second NLCD raster layer with only these three values:
    - 52: "SHRUB/SCRUB"
    - 81: "PASTURELAND"
    - 82: "CROPLAND"
3) Get County and State Shapefiles.
    - We need the State and County FIPS codes and names
4) Rasterize a State + County FIPS code to the same geometry as the NLCD
   raster above:
    - Might need to concatenate these values.
5) Multiply these rasters together:
    - This works because the products create all unique values.
    - Double check this.
6) Associate each product value with its State-County-NLCD combination and
   create a product-combination table.
7) Get a land-value table with price, State-County-NLCD combinations, and
   index values.
8) Join that land value table with the product table.
9) Map the land value table's index to the product raster to create a new
   raster of land value indices.
10) At some point we will need an acre grid of land-values. How will we decide
    which land-values from the 30m grid to assign to acre grid cell?

"""
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import rasterio
import xarray as xr
from functions import rasterize

# Data Path
DP = "/scratch/twillia2/weto_costs/data"

# get the target geometry
nlcd = rasterio.open(os.path.join(DP, "rasters/nlcd_2016_ag.tif"))
geom = nlcd.get_transform()
ny = nlcd.height
nx = nlcd.width
xs = [geom[0] + geom[1] * i  for i in range(nx)]
ys = [geom[3] + geom[-1] * i  for i in range(ny)]
extent = [np.min(xs), np.min(ys), np.max(xs), np.max(ys)]


# The geoid is a combo of state and county fips
if not os.path.exists(os.path.join(DP, "rasters/county_gids.tif")):
    rasterize(src=os.path.join(DP, "shapefiles/USA/tl_2017_us_county.shp"),
              dst=os.path.join(DP, "rasters/county_gids.tif"),
              attribute="GEOID",
              resolution=geom[1],
              cols=nx,
              rows=ny,
              epsg=4326,
              extent=extent,
              overwrite=True)

# We need to combine these values
ag = xr.open_rasterio(os.path.join(DP, "rasters/nlcd_2016_ag.tif"))
county = xr.open_rasterio(os.path.join(DP, "rasters/county_gids.tif"))

# if the product of these two sets of values results in all unique values...
uags = np.array([52, 81, 82])
cdf = gpd.read_file(os.path.join(DP, "shapefiles/USA/tl_2017_us_county.shp"))
ugids = cdf["GEOID"].unique().astype(int)
ufips = cdf["COUNTYFP"].unique().astype(int)

# Are geoid products unique?
product = []
for uag in uags:
    vals = uag * ugids
    vals = list(vals)
    product += vals
assert np.unique(np.array(product)).shape[0] == len(product)  # yes

# Are fips products unique?
product = []
for uag in uags:
    vals = uag * ufips
    vals = list(vals)
    product += vals
assert np.unique(np.array(product)).shape[0] == len(product)  # no


# These are the values we need to associate with
lookup = pd.read_csv(os.path.join(DP, "tables/conus_cbe_lookup.csv"))
lookup.columns = ['index', 'type',' dollar_ac']

# So we need a table with agid (ag + gid) associated values
states = gpd.read_file(os.path.join(DP, "shapefiles/USA/tl_2016_us_state.shp"))
states = states[["STATEFP", "NAME"]]
counties = cdf[["GEOID", "NAME", "NAMELSAD", "STATEFP", "COUNTYFP"]]
reference = pd.merge(counties, states, on="STATEFP")
reference.columns = ['GEOID', 'NAMECTY', 'NAMELSAD', 'STATEFP', 'COUNTYFP',
                     'NAMEST']

# We a column with COUNTY STATE AGTYPE, all caps
capper = lambda x: x["NAMECTY"].upper() + " " + x["NAMEST"].upper()
reference["type"] = reference[["NAMECTY", "NAMEST"]].apply(capper, axis=1)
reference = reference[["GEOID", "type"]]

# Now we need every combination of ag and type (3*n counties)
shrubref = reference.copy()
pastref = reference.copy()
cropref = reference.copy()
shrubref["nlcd"] = 52
pastref["nlcd"] = 81
cropref["nlcd"] = 82
shrubref["legend"] = "SHRUB/SCRUB"
pastref["legend"] = "PASTURELAND"
cropref["legend"] = "CROPLAND"
ref2 = pd.concat([shrubref, cropref, pastref], sort=True)
ref2["type"] = ref2["type"] + " " + ref2["legend"]
ref2["rast_val"] = ref2["GEOID"].astype(int) * ref2["nlcd"].astype(int)
ref2.to_csv(os.path.join(DP, "tables/nlcd_rast_lookup.csv"), index=False)

# Compare these values with the agcounty_product raster
agcounty = xr.open_rasterio(os.path.join(DP, "rasters/agcounty_product.tif"))
val = agcounty[0, 50000, 50000].data
ref2[ref2["rast_val"] == val]
