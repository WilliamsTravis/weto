#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 13:36:28 2020

@author: twillia2
"""

from reVX.utilities.exclusions_converter import ExclusionsConverter
from gdalmethods import Data_Path


DP = Data_Path("/scratch/twillia2/weto/data")
RENT = DP.join("albers/acre/rent_map.tif")

def to_h5(src, dst):
    """Convert the rent GeoTiff to a format that can be used in reV.

    Parameters
    ----------
    src : str
        .
    dst : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    
    # Map to WTK point
    needed = ["wtk gid", "added_costs"]