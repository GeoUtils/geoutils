# -*- coding: utf-8 -*-
"""
georaster.py

Classes to simplify reading geographic raster data and associated 
georeferencing information.

There are two classes available:

    SingleBandRaster    --  For loading datasets with only a single band of 
                            data
    MultiBandRaster     --  For loading datasets that contain multiple bands 
                            of data

Each class works as a wrapper to the GDAL API. A raster dataset can be loaded
into a class without needing to load the actual data as well, which is useful
for querying geo-referencing information without memory overheads.

Both classes provide comprehensive georeferencing information access via
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

Example use for SingleBandRaster with data import:
>>> import georaster
>>> my_image = georaster.SingleBandRaster('myfile.tif')
>>> plt.imshow(my_image.r,extent=my_image.extent)
The data of a SingleBandRaster is made available via my_image.r as a numpy
array.

Example use for SingleBandRaster to get some georeferencing info, without 
also loading data into memory:
>>> import georaster
>>> my_image = georaster.SingleBandRaster('myfile.tif',load_data=False)
>>> print my_image.srs.GetProjParm('central_meridian')

Example use for MultiBandRaster, loading all bands:
>>> import georaster
>>> my_image = georaster.MultiBandRaster('myfile.tif')
>>> plt.imshow(my_image.r)

Example use for MultiBandRaster, loading just a couple of bands:
>>> import georaster
>>> my_image = georaster.MultiBandRaster('myfile.tif',load_data=[1,3])
>>> plt.imshow(my_image.r[:,:,my_image.gdal_band(3)])

For more examples see the the class itself.


Created on Wed Oct 23 12:06:16 2013

@author: Andrew Tedstone (a.j.tedstone@ed.ac.uk)
"""

import numpy as np
from scipy import ndimage
from osgeo import osr, gdal
import subprocess
import xml.etree.ElementTree as etree
import re
import mpl_toolkits.basemap.pyproj as pyproj

# By default, GDAL does not raise exceptions - enable them
# See http://trac.osgeo.org/gdal/wiki/PythonGotchas
gdal.UseExceptions()

class __Raster:
    """
    Attributes:
        ds_file : filepath and name
        ds : GDAL handle to dataset
        extent : extent of raster in order understood by basemap, 
                    [xll,xur,yll,yur], in raster coordinates
        srs : OSR SpatialReference object
        proj : pyproj conversion object raster coordinates<->lat/lon 

    """      
    # Filepath and name
    ds_file = None
    # GDAL handle to the dataset
    ds = None
    # Extent of raster in order understood by Basemap
    extent = None
    # SRS
    srs = None
    # pyproj Projection
    proj = None
  


    def __del__(self):
        """ Close the gdal link to the dataset on object destruction """
        self.ds = None



    def _load_ds(self,ds_filename):
        """ Load link to data file and set up georeferencing """
        self.ds_file = ds_filename
        self.ds = gdal.Open(ds_filename)
        trans = self.ds.GetGeoTransform()
        self.extent = (trans[0], trans[0] + self.ds.RasterXSize*trans[1], trans[3] + self.ds.RasterYSize*trans[5], trans[3])
        self.srs = osr.SpatialReference()
        self.srs.ImportFromWkt(self.ds.GetProjection())
        self.proj = pyproj.Proj(self.srs.ExportToProj4())


 
    def get_extent_latlon(self):
        """ Return raster extent in lat/lon, (xll,xur,yll,yur) """
        left,bottom = self.proj(self.extent[0],self.extent[2],inverse=True)
        right,top = self.proj(self.extent[1],self.extent[3],inverse=True)
        return (left,right,bottom,top)



    def read_single_band(self,band=1):
        """ Return np array of specified band, defaults to band 1. """
        return self.ds.GetRasterBand(band).ReadAsArray()  



    def read_single_band_subset(self,bounds,latlon=False,band=1):
        """ Return a subset area of the specified band of the dataset.

        Supply coordinates in native system or lat/lon.

        Parameters:
            bounds : tuple (xstart,xend,ystart,yend) (same format as extent)
            latlon : boolean, default False. Set as True if bounds in lat/lon.
            band : int, number of band to read. Default 1.

        Returns:
            np.array of subset area

        Map --> px conversion from:
        http://gis.stackexchange.com/a/46898

        """
        
        left,right,bottom,top = bounds

        # Convert coordinates to map system if provided in lat/lon
        if latlon == True:
            left,bottom = self.proj(left,bottom)
            right,top = self.proj(right,top)

        # Start conversion to pixel coordinate space.
        # Unlike the bounds tuple, which specifies bottom left and top right
        # coordinates, here we need top left and bottom right for the numpy
        # readAsArray implementation.
        trans = self.ds.GetGeoTransform()
        # Coordinates of start point in pixel coordinates
        xpx1 = int((left - trans[0]) / trans[1])
        ypx1 = int((top - trans[3]) / trans[5]) 
        # Coordinates of finish point in pixel coordinates
        xpx2 = int((right - trans[0]) / trans[1])
        ypx2 = int((bottom - trans[3]) / trans[5])
        # Resulting pixel offsets
        x_offset = xpx2 - xpx1
        y_offset = ypx2 - ypx1
        # In special case of being called to read a single point, offset 1 px
        if x_offset == 0: x_offset = 1
        if y_offset == 0: y_offset = 1

        # Read array and return
        arr = self.ds.GetRasterBand(band).ReadAsArray(xpx1,ypx1,
                                                      x_offset,y_offset)
        return arr



    def value_at_coords(self,x,y,latlon=False,band=None,system=None):
        """ Extract the pixel value(s) at the specified coordinates.
        
        Extract pixel value of each band in dataset at the specified 
        coordinates. Alternatively, if band is specified, return only that
        band's pixel value.
                
        Parameters:
            x : float, x coordinate.
            y : float, y coordinate.
            latlon : boolean, True if coordinates in WGS84, false if in 
                     native system of the raster.
            band : the band number to extract from.
            system : DEPRECATED but maintained for backward compatibility.

        Returns:
            if band specified, float of extracted pixel value.

            if band not specified, depending on number of bands in dataset:
                more than 1 band : dict of GDAL band number:float pixel value.
                just 1 band : float of extracted pixel value.
        
        """

        def format_value(value):
            """ Check if valid value has been extracted """
            if type(value) == np.ndarray:
                value = value[0,0]
            else:
                value = None
            return value

        # Check the deprecated georeferencing system parameter.
        if system == 'wgs84':
            latlon = True

        # Create instance of bounds to read; read_single_band_subset will
        # increment these by 1px after completing the necessary georeferencing
        # operations.
        bounds = [x,x,y,y]

        # N.b. Values are returned in array from read_single_band_subset so we 
        # index in to retrieve the single coordinate value.

        # Get values for all bands
        if band == None:

            # Deal with SingleBandRaster case
            if self.ds.RasterCount == 1:
                value = format_value(self.read_single_band_subset(bounds,
                                                            latlon=latlon))
            
            # Deal with MultiBandRaster case
            else:    
                value = {}
                for b in range(1,self.ds.RasterCount+1):
                    val = format_value(self.read_single_band_subset(bounds,
                                                    latlon=latlon,band=b))
                    # Store according to GDAL band numbers
                    value[b] = val

        # Or just for specified band in MultiBandRaster case                
        elif isinstance(band,int):
            value = format_value(self.read_single_band_subset(bounds,
                                                    latlon=latlon,band=band))
        else:
            raise ValueError('The value provided for band was not int or None.')

        return value

    

    def coordinates(self):
        """ Calculate x and y coordinates for every pixel. 

        Parameters:
            None.
        Returns:
            (x,y) : tuple, containing 1-d numpy array for each of x and y.
            
        """
        trans = self.ds.GetGeoTransform()
        if self.ds.RasterCount > 1:
            shape = self.r.shape[0:2]
        else:
            shape = self.r.shape
        x = np.array(np.linspace(trans[0],(trans[0]+(self.ds.RasterXSize*trans[1])),self.ds.RasterXSize).tolist() * self.ds.RasterYSize).reshape(shape)
        y = np.array(np.linspace(trans[3],(trans[3]+(self.ds.RasterYSize*trans[5])),self.ds.RasterYSize).tolist() * self.ds.RasterXSize).reshape(shape[::-1]).T
        return (x,y)



    def get_utm_zone(self):
        """ 
        Return UTM zone of raster from GDAL Projection information. 

        This function used to be more complex but is now a wrapper to an OGR
        Spatial Reference call. It remains maintained for backwards 
        compatibility with dependent scripts.

        Returns:
            str : zone

        """
        return self.srs.GetUTMZone()
        


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

    Attributes:
        ds_file : filepath and name
        ds : GDAL handle to dataset
        extent : extent of raster in order understood by basemap, 
                    [xll,xur,yll,yur], in raster coordinates
        srs : OSR SpatialReference object
        proj : pyproj conversion object raster coordinates<->lat/lon 
        r : numpy band of array data

    Example:
    >>> georaster.SingleBandRaster('myfile.tif',load_data=True|False)

    """
     
    # Numpy array of band data
    r = None
     
    def __init__(self,ds_filename,load_data=True):
        """ Construct object with raster from a single band dataset. 
        
        Parameters:
            ds_filename : filename of the dataset to import
            load_data : boolean, default True, to import the data into obj.r.
            
        """
        self._load_ds(ds_filename)     
        
        if load_data == True:
            self.r = self.read_single_band(1)

    
        
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
                


    def smooth(self,px=3):
        """ Apply gaussian filter of specified window px.
        
        Parameters:
            px : window size in pixels to supply to ndimage.gaussian_filter
        Returns:
            self.smoothed : smoothed version of self.r.
            
        """
        self.smoothed = ndimage.gaussian_filter(self.r,px)
        



class MultiBandRaster(__Raster):
    """ A geographic raster dataset with multiple bands of data.
    
    Initialise with the file path to a single band raster dataset, of type
    understood by the GDAL library. Datasets in UTM are preferred.

    Examples:
    [1] Load all bands of a raster:
    >>> georaster.MultiBandRaster("myfile.tif")

    [2] Load just bands 1 and 3 of a raster:
    >>> georaster.MultiBandRaster("myfile.tif",load_data=[1,3])

    [3] Don't load the data, just use the class for the georeferencing API:
    >>> georaster.MultiBandRaster("myfile.tif",load_data=False) 

    Attributes:
        r :         The raster data (if loaded). For the example cases above:

                    [1] - a np array of [rows,cols,bands]. Standard numpy 
                          slicing is used, ie. the array is zero-referenced.
                          E.g., extract band 2, which is the second band 
                          loaded:
                          >>> myRaster.r[:,:,1]

                          This can be simplified with gdal_band():
                          E.g., extract band 2:
                          >>> myRaster.r[:,:,myRaster.gdal_band(2)]

                    [2] - The same as [1]. The helper function is particularly 
                          useful in simplifying lookup of bands, e.g.:
                          >>> myRaster.r[:,:,myRaster.gdal_band(3)]
                          Rather than the less obvious:
                          >>> myRaster.r[:,:,1]
                          Which corresponds to the actual numpy location of 
                          that band.

                    [3] - r is set as None. No data can be accessed.


        bands :     list of GDAL band numbers which have been loaded, in the 
                    order corresponding to the order stored in 
                    r[rows,cols,bands].  

        ds_file :   filepath and name
        ds :        GDAL handle to dataset
        extent :    extent of raster in order understood by basemap, 
                    [xll,xur,yll,yur], in raster coordinates
        srs :       OSR SpatialReference object
        proj :      pyproj conversion object raster coordinates<->lat/lon 
        
    """

    # Either a numpy array if just one band, or a dict of numpy arrays if 
    # multiple, key is band number
    r = None

    # List of GDAL band numbers which have been loaded.
    bands = None


    def __init__(self,ds_filename,load_data=True):
        """ Load a multi-band raster.

        Parameters:
            ds_filename : filename of dataset to load
            load_data : True, False, or tuple of raster bands. If tuple,
                MultiBandRaster.r will be a numpy array [y,x,b], where bands
                are indexed from 0 to n in the order specified in the tuple.

        """
        self._load_ds(ds_filename)

        if load_data == True:
            self.r = np.zeros((self.ds.RasterYSize,self.ds.RasterXSize,
                               self.ds.RasterCount))
            for b in range(0,self.ds.RasterCount):
                self.r[:,:,b] = self.read_single_band(band=b+1)
            self.bands = range(1,self.ds.RasterCount+1)

        elif load_data == False:
            bands = None
            return

        elif isinstance(load_data,(tuple,list)) == True:
            self.r = np.zeros((self.ds.RasterYSize,self.ds.RasterXSize,
                               len(load_data)))
            k = 0
            for b in load_data:
                self.r[:,:,k] = self.read_single_band(band=b)
                k += 1
            self.bands = load_data

        else:
            raise ValueError('load_data was not understood (should be one \
             of True,False or tuple)')



    def gdal_band(b):
        """ Return numpy array location index for given GDAL band number. 

        Parameters:
            b : int, value of GDAL band number to lookup

        Returns:
            int, index location of band in self.r

        Example:
        >>> giveMeMyBand2 = myRaster.r[:,:,myRaster.gdal_band(2)]

        """

        # Check that more than 1 band has been loaded into memory.
        if self.bands == None:
            raise AttributeError('No data have been loaded.')
        if len(self.bands) == 1:
            raise AttributeError('Only 1 band of data has been loaded.')

        if isinstance(b,int):
            return self.bands.index(b)
        else:
            raise ValueError('B is must be an integer.') 







def simple_write_geotiff(outfile,raster,geoTransform,wkt=None,proj4=None):
    """ Save a GeoTIFF.
    
    Inputs:
        outfile : filename to save image to
        raster : nbands x r x c
        geoTransform : tuple (top left x, w-e cell size, 0, top left y, 0, n-s cell size (-ve))
        One of proj4 or wkt :
            proj4 : a proj4 string
            wkt : a WKT projection string

    Outputs:
        A GeoTiff named outfile.
    
    Based on http://adventuresindevelopment.blogspot.com/2008/12/python-gdal-adding-geotiff-meta-data.html
    and http://www.gdal.org/gdal_tutorial.html
    """

    # Georeferencing sanity checks
    if wkt <> None and proj4 <> None:
        raise 'InputError: Both wkt and proj4 specified. Only specify one.'
    if wkt == None and proj4 == None:
        raise 'InputError: One of wkt or proj4 need to be specified.'

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
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(outfile, xdim, ydim, nbands, gdal.GDT_Float32)
    # Top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    dst_ds.SetGeoTransform(geoTransform)
      
    # Set the reference info 
    srs = osr.SpatialReference()
    if wkt <> None:
        dst_ds.SetProjection(wkt)
    elif proj4 <> None:
        srs.ImportFromProj4(proj4)
        dst_ds.SetProjection( srs.ExportToWkt() )
    
    # Write the band(s)
    if nbands > 1:
        for band in range(1,nbands+1):
            dst_ds.GetRasterBand(band).WriteArray(raster[band-1]) 
    else:
        dst_ds.GetRasterBand(1).WriteArray(raster)

    # Close data set
    dst_ds = None
    return True 







def write_geotiff(raster,output_file,geo,proj4=False,wkt=False,mask=None):
    """ Save a GeoTIFF. DEPRECATED - DO NOT USE FOR NEW APPLICATIONS
    
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


        
        
