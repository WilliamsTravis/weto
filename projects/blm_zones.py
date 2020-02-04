#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rasterize BLM zone and state land lookup codes.

1) Get a US county shapefile
2) Either use a pdf table converter or find a new version of the 
2) Join the code value to each zone from our lookup.
3) Use the geometry of the WGS NLCD codes GeoTiff to rasterize
4) Reproject into Albers Equal Area Conic (EPSG 102008)
5) I think if we get the 30m raster, the hard part will be done and we can
   resample to whatever resolution later.


To Do:
    It looks like these cover all the land in the BLM territory, so we'll have
    clip out just the fed land portion first.
    
Created on Mon Feb  3 09:30:06 2020

@author: twillia2
"""
import camelot
import geopandas as gpd
import os
import pandas as pd
import rasterio
import requests
from osgeo import gdal
from weto.functions import Data_Path, rasterize


# Data Path
dp = Data_Path("/scratch/twillia2/weto/data")

# Lookup table
lookup = pd.read_csv(dp.join("tables/conus_cbe_lookup.csv"))
lookup.columns = ['code', 'type',' dollar_ac']
zone_lu = lookup[lookup["type"].str.contains("BLM Zone")]
zone_lu["zone"] = zone_lu["type"].apply(lambda x: int(x[-2:]))


# Tabula py!
pdf_path = dp.join("tables/blm_rent_zones.pdf")
if not os.path.exists(pdf_path):
    url = ("https://www.blm.gov/download/file/fid/18524/IM%202017-096%20" +
           "Attachment%205%20Adjusted%202012%20NASS%20Census%20Per%20Acre%" +
           "20LB%20Values%20and%20Rent%20Schedule%20Zones.pdf")
    r = requests.get(url)
    with open(pdf_path, 'wb') as file:
        file.write(r.content)

# Get the table from each page - this mostly works
pdfs = camelot.read_pdf(pdf_path, pages="1-end")
dfs = [pdfs[i].df for i in range(pdfs.n)]
zone_df = pd.concat(dfs, sort=False).reset_index(drop=True)
headers = zone_df.loc[0]
zone_df = zone_df.iloc[1:]
zone_df = zone_df.reset_index(drop=True)
zone_df.columns = ["state", "county", "pfactor", "value", "zone"]

# To join we need to clean up the county names
fixit = lambda x: x.replace("*", "").replace("\n", "")
zone_df = zone_df[zone_df["state"] != "Alaska"]
zone_df["county"] = zone_df["county"].apply(fixit)

# We need a state county field
counties = gpd.read_file(dp.join("shapefiles/USA/tl_2017_us_county.shp"))
cols = list(counties.columns[:6] )+ ["geometry"]
counties = counties[cols]
counties.columns = ['STATEFP', 'COUNTYFP', 'COUNTYNS', 'GEOID', 'NAME',
                    'NAMELSAD', 'geometry']
states = gpd.read_file(dp.join("shapefiles/USA/tl_2017_us_state.shp"))
states = states[["NAME", "STATEFP", "geometry"]]
counties = counties[["NAME", "STATEFP", "geometry"]]
zone_df["stateco"] = zone_df["county"] + ", " + zone_df["state"]
counties = counties.merge(states[["NAME", "STATEFP"]], on="STATEFP")
counties["stateco"] = counties["NAME_x"] + ", " + counties["NAME_y"]
counties = counties.drop(["NAME_x", "NAME_y"], axis=1)

# Now merge everything together
zone_df = zone_df[zone_df["zone"] != ""]
zone_df["zone"] = zone_df["zone"].astype(int)
zone_df = zone_df.merge(zone_lu, on="zone", how="left")
counties = counties.merge(zone_df, on="stateco", how="left")
counties = counties[['STATEFP', 'geometry', 'stateco', 'state', 'county',
                     'pfactor', 'value', 'zone', 'code', 'type',
                     ' dollar_ac']]

# Now clip by fedlands - do this with ogr it takes too long
fedland = gpd.read_file(dp.join("shapefiles/BLM/conus_fedland_blm_county.shp"))
fedland = fedland[["NAMELSAD", "geometry"]]
blm_counties = gpd.overlay(counties, fedland)

# Save to new file
blm_counties.to_file(dp.join("shapefiles/USA/blm_2017_us_county.shp"))

# get the target geometry
nlcd = rasterio.open(dp.join("rasters/nlcd_2016_ag.tif"))
transform = nlcd.get_transform()
height = nlcd.height
width = nlcd.width

# Rasterize the code field to the finest resolution we'll use
rasterize(src=dp.join("shapefiles/USA/blm_2017_us_county.shp"),
          dst=dp.join("rasters/blm_codes.tif"),
          attribute="code",
          transform=transform,
          height=height,
          width=width,
          epsg=4326,
          dtype=gdal.GDT_Float32,
          overwrite=True)
