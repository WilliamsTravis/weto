# -*- coding: utf-8 -*-
"""
Create coverage statistics for costs.

Created on Tue Mar 24 10:19:58 2020

@author: twillia2
"""

import pandas as pd
from coverage_codes import get_counts, coverage
from gdalmethods import Data_Path


# DP = Data_Path("~/data/weto/rent_map")
DP = Data_Path("/scratch/twillia2/weto/data")


def fixit(x):
    x = x.replace("*", "").replace("\n", "").replace("$", "").replace(",", "")
    if x == " NULL " or  x == "  NULL  ":
        x = 0
    return x


def codes():
    """Collect the cost codes from each category into a dictionary."""

    # Code to cost lookup table
    lookup = pd.read_csv(DP.join("tables/conus_cbe_lookup.csv"))
    lookup.columns = ['code', 'type', 'dollar_ac']
    lookup["dollar_ac"] = lookup["dollar_ac"].apply(fixit)
    lookup["dollar_ac"] = lookup["dollar_ac"].astype(float)

    # Now, change the codes to only include those with values
    lookup = lookup[lookup["dollar_ac"] > 0.0]

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


if __name__ == "__main__":
    codes = codes()
    counts = get_counts(codes)
    covdf = coverage(counts)
