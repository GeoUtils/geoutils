# -*- coding: utf-8 -*-
"""
georaster.py

Created on Wed Oct 23 12:06:16 2013

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
from scipy import ndimage
from osgeo import osr, gdal
import subprocess
import xml.etree.ElementTree as etree
import re

class __Raster:
    
    # Filepath and name
    ds_file = None
    # GDAL handle to the dataset
    ds = None
    # Extent of raster in order understood by Basemap
    extent = None
    
    def __del__(self):
        """ Close the gdal link to the dataset on object destruction """
        self.ds = None


    def read_single_band(self,band):
        """ Return np array of specified band. """
        return self.ds.GetRasterBand(band).ReadAsArray()  


    def value_at_coords(self,x,y,band=None):
        """ Extract the pixel value(s) at the specified coordinates.
        
        Extract pixel value of each band in dataset at the specified 
        coordinates. Alternatively, if band is specified, return only that
        band's pixel value.
                
        Parameters:
            x : x coordinate in format of target dataset.
            y : y coordinate in format of target dataset.
            band : the band number to extract from.
        Returns:
            if band specified, float of extracted pixel value.
            if band not specified, dict of band number:float pixel value.
        
        """
        cmd = "gdallocationinfo " + self.ds_file + " -xml -geoloc " + str(x) + " " + str(y)  
        pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        stdout,stderr = pid.communicate()
        xmlr = etree.fromstring(stdout) 
        
        vals = {}
        for r in xmlr.findall('BandReport'):
            b = r.attrib['band']
            v = float(r.find('Value').text)
            if b == band:
                return v
            else:
                vals[b] = v
            
        return vals

    
    def coordinates(self):
        """ Calculate x and y coordinates for every pixel. 

        Parameters:
            None.
        Returns:
            (x,y) : tuple, containing 1-d numpy array for each of x and y.
            
        """
        trans = self.ds.GetGeoTransform()
        x = np.array(np.linspace(trans[0],(trans[0]+(self.ds.RasterXSize*trans[1])),self.ds.RasterXSize).tolist() * self.ds.RasterYSize).reshape(self.r.shape)
        y = np.array(np.linspace(trans[3],(trans[3]+(self.ds.RasterYSize*trans[5])),self.ds.RasterYSize).tolist() * self.ds.RasterXSize).reshape(self.r.shape[::-1]).T
        return (x,y)


    def get_utm_zone(self):
        """ 
        Return UTM zone of raster from GDAL Projection information. 

        Returns:
            str : zone

        """
        value = re.findall('UTM zone ([0-9]{1,2}N|S)',self.ds.GetProjection())
        if len(value) <> 1:
            raise 'UTM zone not found'
        return value[0]


    def get_pixel_size(self):
        """ 
        Return pixel size of loaded raster.

        Returns:
            floats: xres, yres

         """
        geotransform = self.ds.GetGeoTransform()
        xres = geotransform[1]
        yres = geotransform[5]
        return xres, yres






class SingleBandRaster(__Raster):
    """ A geographic raster dataset with one band of data.
    
    Initialise with the file path to a single band raster dataset, of type
    understood by the GDAL library. Datasets in UTM are preferred.
    
    """
     
    # Numpy array of band data
    r = None
     
    def __init__(self,ds_filename,load_data=True):
        """ Construct object with raster from a single band dataset. 
        
        Parameters:
            ds_filename : filename of the dataset to import
            load_data : boolean, default True, to import the data into obj.r.
        Returns:
            Sets self.ds, self.extent, self.r.
            
        """
        self.ds_file = ds_filename
        self.ds = gdal.Open(ds_filename)
        trans = self.ds.GetGeoTransform()
        self.extent = (trans[0], trans[0] + self.ds.RasterXSize*trans[1], trans[3] + self.ds.RasterYSize*trans[5], trans[3])
        
        if load_data == True:
            self.r = self.read_single_band(1)

    
        
    def find_value_at_coords(self,x,y):
        """ Extract the pixel value at the specified coordinates.
        
        Parameters:
            x : x coordinate in format of target dataset.
            y : y coordinate in format of target dataset.
        Returns:
            float of extracted pixel value.
        
        """
        return self.value_at_coords(x,y,band=1)
             
        
    def transect(self):
        """ Not yet started.
        
        """
        
    def smooth(self,px=3):
        """ Apply gaussian filter of specified window px.
        
        Parameters:
            px : window size in pixels to supply to ndimage.gaussian_filter
        Returns:
            self.smoothed : smoothed version of self.r.
            
        """
        self.smoothed = ndimage.gaussian_filter(self.r,px)
        


class MultiBandRaster(__Raster):
    """ Placeholder """



def write_geotiff(raster,output_file,geo,proj4=False,wkt=False,mask=None):
    """ Save a GeoTIFF.
    
    Inputs:
        raster - nbands x r x c
        output_file - filename to save image to
        geo - georeferencing information as dictionary. x (left coord), 
        y (top coord), xcellsize, ycellsize, datum, utmzone (if raster is projected in UTM)
        proj4 - a proj4 string, optional
        wkt - a WKT projection string, optional
        Only provide one of proj4 or wkt!

    Outputs:
        A GeoTiff named according to output_file.
    
    Based on http://adventuresindevelopment.blogspot.com/2008/12/python-gdal-adding-geotiff-meta-data.html
    and http://www.gdal.org/gdal_tutorial.html
    """   
    # Check if the image is multi-band or not. 
    if raster.shape.__len__() == 3:
        nbands = raster.shape[0]    
        ydim = raster.shape[1]
        xdim = raster.shape[2]
    elif raster.shape.__len__() == 2:
        nbands = 1
        ydim = raster.shape[0]
        xdim = raster.shape[1]
         
    # Setup geotiff file.
    gdal.SetConfigOption('GDAL_TIFF_INTERNAL_MASK', 'YES')
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_file, xdim, ydim, nbands, gdal.GDT_Float32)
    dst_ds.CreateMaskBand(gdal.GMF_PER_DATASET)
    # Top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    dst_ds.SetGeoTransform( [ geo['x'], geo['xcellsize'], 0, geo['y'], 0, geo['ycellsize'] ] )
      
    # Set the reference info 
    srs = osr.SpatialReference()
    if wkt == False:
        srs.ImportFromProj4(proj4)

        if geo.has_key('utmzone') == True:
            srs.SetUTM(geo['utmzone'],1)
        srs.SetWellKnownGeogCS(geo['datum'])
        dst_ds.SetProjection( srs.ExportToWkt() )
    else:
        dst_ds.SetProjection(wkt)
    
    # Write the band(s)
    if nbands > 1:
        for band in range(1,nbands+1):
            dst_ds.GetRasterBand(band).WriteArray(raster[band-1]) 
            if mask <> None:
                dst_ds.GetRasterBand(band).GetMaskBand().WriteArray(mask)
    else:
        dst_ds.GetRasterBand(1).WriteArray(raster)
        if mask <> None:
            dst_ds.GetRasterBand(1).GetMaskBand().WriteArray(mask)

    # Close data set
    dst_ds = None
    return True 
        
        
        
