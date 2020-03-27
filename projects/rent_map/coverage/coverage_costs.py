# -*- coding: utf-8 -*-
"""
Create coverage statistics for costs.

Created on Tue Mar 24 10:19:58 2020

@author: twillia2
"""

from coverage_codes import get_counts, get_coverage


if __name__ == "__main__":
    counts = get_counts(cost_coverage=True)
    covdf = get_coverage(counts, file_path="coverage_cost.csv")
