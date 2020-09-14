#!/usr/bin/env python

from distutils.core import setup

#To prepare a new release
#python setup.py sdist upload

setup(name='geoutils',
    version='1.0',
    description='Libraries and command-line utilities for geospatial data processing/analysis',
    author='Amaury Dehecq',
    author_email='amaury.dehecq@jpl.nasa.gov',
    license='MIT',
      url='https://bitbucket.org/atedstone/geoutils.git',
    packages=['geoutils'],
    #long_description=open('README.md').read(),
    #install_requires=['gdal','numpy','scipy','matplotlib'],
    #Note: this will write to /usr/local/bin
    scripts=['geoutils/crop2common.py','geoutils/crop2image.py','geoutils/dem_coregistration.py','geoutils/plot_geo_image.py'])

