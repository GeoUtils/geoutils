# -*- coding: utf-8 -*-
"""
demraster.py

Classes to simplify reading DEM rasters and associated 
georeferencing information.
Similar to georaster classes with extra functions for DEM display/processing

The class provides comprehensive georeferencing information access via
the following attributes:

    class.ds    :   the GDAL handle to the dataset, which provides access to 
                    all GDAL functions, e.g. GetProjection, GetGeoTransform.
                        More information on API:
                        http://www.gdal.org/classGDALDataset.html

    class.srs   :   an OGR Spatial Reference object representation of the 
                    dataset.
                        More information on API:
                        http://www.gdal.org/classOGRSpatialReference.html

    class.proj  :   a pyproj coordinate conversion function between the 
                    dataset coordinate system and lat/lon.

    class.extent :  tuple of the corners of the dataset in native coordinate
                    system, as (left,right,bottom,top).

Additionally, georeferencing information which requires calculation 
can be accessed via a number of functions available in the class.


Raster access:
    
    class.r : numpy array, m*n for SingleBandRaster


Created on Mon Apr 6

@author: Amaury Dehecq (amaury.dehecq@univ-savoie.fr)
"""

import numpy as np
from scipy import ndimage, signal
from scipy.stats import nanmean
from osgeo import osr, gdal
import subprocess
import xml.etree.ElementTree as etree
import re
import mpl_toolkits.basemap.pyproj as pyproj

from georaster import __Raster
import geometry as geo

# By default, GDAL does not raise exceptions - enable them
# See http://trac.osgeo.org/gdal/wiki/PythonGotchas
gdal.UseExceptions()


class DEMRaster(__Raster):
    """ A geographic raster dataset with one band of data.
    
    Initialise with the file path to a single band raster dataset, of type
    understood by the GDAL library. Datasets in UTM are preferred.

    Attributes:
        ds_file : filepath and name
        ds : GDAL handle to dataset
        extent : extent of raster in order understood by basemap, 
                    [xll,xur,yll,yur], in raster coordinates
        srs : OSR SpatialReference object
        proj : pyproj conversion object raster coordinates<->lat/lon 
        r : numpy band of array data

    Example:
    >>> demraster.DEMRaster('myfile.tif',load_data=True|False)

    """
     
    # Numpy array of band data
    r = None
     
    def __init__(self,ds_filename,load_data=True,latlon=True,band=1):
        """ Construct object with raster from a single band dataset. 
        
        Parameters:
            ds_filename : filename of the dataset to import
            load_data : - True, to import the data into obj.r. 
                        - False, to not load any data.
                        - tuple (left, right, bottom, top) to load subset; 
                          obj.extent will be set to reflect subset area.
            latlon : default True. Only used if load_data=tuple. Set as False
                     if tuple is projected coordinates, True if WGS84.
            band : default 1. Specify GDAL band number to load. If you want to
                   load multiple bands at once use MultiBandRaster instead.
            
        """
        self._load_ds(ds_filename)     
        
        # Load entire image
        if load_data == True:
            self.r = self.read_single_band(band)

        # Or load just a subset region
        elif isinstance(load_data,tuple) or isinstance(load_data,list):
            if len(load_data) == 4:
                (self.r,self.extent) = self.read_single_band_subset(load_data,
                                        latlon=latlon,extent=True,band=band)

        elif load_data == False:
            return

        else:
            print 'Warning : load_data argument not understood. No data loaded.'

   
        
    def find_value_at_coords(self,x,y,**kwargs):
        """ (DEPRECATED) Extract the pixel value at the specified coordinates.
        
        This function is maintained for backward compatibility. New code 
        should call self.value_at_coords() instead.

        ----------------------------------------------------------------------

        Parameters:
            x : x coordinate in format of target dataset.
            y : y coordinate in format of target dataset.
        Returns:
            float of extracted pixel value.
        
        """
        return self.value_at_coords(x,y,**kwargs)
                


    def compute_slope(self):
        '''
        Calculate the slope and aspect of a DEM
        Outputs :
        - slope : matrix of same size as input DEM, slope in radians
        - aspect : matrix of same size as input DEM, aspect in radians
        '''
        
        #compute distance in meters with right/lower neighbour
        if self.srs.IsProjected()==0:
            lon, lat = self.coordinates()
            lonr = np.roll(lon,1,1)
            latl = np.roll(lat,1,0)
            distx = geo.dist_ortho(lon,lat,lonr,lat)
            disty = geo.dist_ortho(lon,lat,lon,latl)
        else:
            print "Need implementation for projected systems"
            sys.exit(1)

        #Compute z gradients
        f1 = np.array([[-1,0,1],[-2,0,2],[-1,0,1]])
        f2 = f1.transpose()
        g1 = signal.convolve(self.r,f1,mode='same')
        g2 = signal.convolve(self.r,f2,mode='same')

        #compute slope
        mpm = np.sqrt((g1/distx)**2 + (g2/disty)**2)/8
        slope = np.arctan(mpm)
        slope_deg = (slope*180)/(np.pi)
        aspect = np.arctan2(g2,g1)

        return slope, aspect


