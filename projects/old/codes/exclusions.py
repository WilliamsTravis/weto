#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
We need a mask for just full exclusions. 1s for full exclusions 0s for
everything else.

Created on Fri Feb 14 09:18:28 2020

@author: twillia2
"""

import geopandas as gpd
import os
import pandas as pd
import rasterio
from osgeo import gdal
from weto.functions import Data_Path, rasterize

# Data Path
dp = Data_Path("/scratch/twillia2/weto/data")
