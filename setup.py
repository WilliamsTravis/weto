# -*- coding: utf-8 -*-
"""
Install methods for calculating cost surfaces for the LandBOSSE wind production
model. 

Created on Thu Jan 30 08:50:16 2020

@author: twillia2
"""

from setuptools import setup

setup(
    name="weto",
    version="0.0.1",
    packages=["weto"],
    description=("Methods for calculating cost surfaces for the LandBOSSE " +
                 "wind production model"),
    author="Travis Williams",
    author_email="travis.williams@nrel.gov",
    install_requires=["dask", "dask-jobqueue", "descartes", "distributed",
                      "geopandas", "h5py", "rasterio", "scipy", "tqdm",
		      "xarray"]
    )
