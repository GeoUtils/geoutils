# Libraries and command-line utilities for geospatial data processing/analysis in Python #

**DEPRECATION NOTICE:** As of 19 February 2021, this repository is marked deprecated and read-only. No further maintenance will be undertaken. Our efforts are now focused on https://www.github.com/GlacioHack/GeoUtils , which incorporates a lot of the functionality of GeoUtils but is built on top of rasterio and geopandas rather than GDAL/OGR directly. Check it out!

---

This package offers Python classes and functions as well as command line tools to work with both geospatial raster and vector datasets. It is built upon the Geospatial Data Abstraction Library (GDAL) bindings, and so in a single command can import any geo-referenced dataset that is understood by GDAL/OGR, complete with all geo-referencing information and various helper functions.

# Full Documentation #

http://geoutils.readthedocs.io


# Installation #

* Summary of set up

`pip install -e .` or `python setup.py install`

* Dependencies

Regular Python packages: numpy, scipy, matplotlib... + GDAL

## Get in touch

* Report bugs, suggest features or view the code on GitHub.
