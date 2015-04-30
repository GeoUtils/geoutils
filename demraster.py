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
import matplotlib.pyplot as pl

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

        #Set nodata values to Nan
        band=self.ds.GetRasterBand(1)
        nodata=band.GetNoDataValue()
        if nodata != None:
            self.r[self.r==nodata] = np.nan


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



    def shaded_relief(self,cmap=pl.get_cmap('Greys'),alpha=1,azdeg=100,altdeg=65,aspect='default'):
        """
        Function to plot a shaded relief of the DEM
        cmap : object of class Colormap
        alpha : f, 0 to 1, transparency
        azdeg : f, azimuth (measured clockwise from south) of the light source in degress for the shaded relief
        altdeg : f, altitude (measured up from the plane of the surface) of the light source in degrees
        aspect : f, imshow window aspect, default is set to have orthonormal axis (meters not degree)
        """

        #set aspect so that lat/lon spacing are equal
        if aspect=='default':
            lat0 = self.extent[2]
            aspect = 1/np.cos(np.abs(lat0)*np.pi/180)

        
        #compute shaded image
        x, y = np.gradient(self.r)  
        slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))  
        aspect = np.arctan2(-x, y)  
        azrad = azdeg*np.pi / 180.  
        altituderad = altdeg*np.pi / 180.  
        shaded = np.sin(altituderad) * np.sin(slope) + np.cos(altituderad) * np.cos(slope)* np.cos(azrad - aspect)  
        rgb=255*(shaded + 1)/2  

        #plot
        pl.imshow(rgb,extent=self.extent,interpolation='bilinear',cmap=cmap,alpha=alpha,aspect=aspect)

        #other option (more time/memory consuming)
        # from matplotlib.colors import LightSource
        # ls = LightSource(azdeg=azdeg,altdeg=altdeg)
        # rgb = ls.shade(self.r,cmap=cm)  
        # pl.imshow(rgb,extent=self.extent,interpolation='bilinear',aspect=aspect,alpha=alpha)
    
            

    def plot_contours(self,levels='default',c='k',alpha=1,aspect='default'):
        """
        Function to plot a shaded relief of the DEM
        cmap : object of class Colormap
        alpha : f, 0 to 1, transparency
        azdeg : f, azimuth (measured clockwise from south) of the light source in degress for the shaded relief
        altdeg : f, altitude (measured up from the plane of the surface) of the light source in degrees
        aspect : f, imshow window aspect, default is set to have orthonormal axis (meters not degree)
        """

        #set aspect so that lat/lon spacing are equal
        if aspect=='default':
            lat0 = self.extent[2]
            aspect = 1/np.cos(np.abs(lat0)*np.pi/180)

        #default level values
        if levels=='default':
            min = np.percentile(self.r,10)
            max = np.nanmax(self.r)
            levels = np.linspace(min,max,10)
            levels = levels//100*100

        #plot
        CS = pl.contour(self.r,levels,extent=self.extent,colors=c,alpha=alpha,aspect=aspect)
        pl.clabel(CS, inline=0, fontsize=10,fmt="%.0f")
        
            

