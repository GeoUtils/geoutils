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
from scipy.ndimage.filters import gaussian_filter
from matplotlib.colors import LightSource

from georaster import __Raster
import geometry as geo
import EGM96

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
     
    def __init__(self,ds_filename,load_data=True,latlon=True,band=1,ref='EGM96'):
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
            ref : EGM96 or WGS84, SRTM data are referenced to the geoid EGM96 while CS2 and SAR data are referenced to the ellipsoid WGS84

        """
        self._load_ds(ds_filename)     
        
        # Load entire image
        if load_data == True:
            self.r = self.read_single_band(band)

        # Or load just a subset region
        elif isinstance(load_data,tuple) or isinstance(load_data,list):
            if len(load_data) == 4:
                (self.r,self.extent) = self.read_single_band_subset(load_data,
                                                                    latlon=latlon,extent=True,band=band,update_info=True)
                
        elif load_data == False:
            return

        else:
            print 'Warning : load_data argument not understood. No data loaded.'

        # Convert to float32 for future computations
        self.r = np.float32(self.r)
        
        #Set nodata values to Nan
        band=self.ds.GetRasterBand(1)
        nodata=band.GetNoDataValue()
        if nodata != None:
            self.r[self.r==nodata] = np.nan

        #Correction to the reference ellipsoid
        self.ref = ref
        if ref=='WGS84':
            if self.srs.IsProjected():
                x, y = self.coordinates()
                lons, lats = self.proj(x,y,inverse=True)
            else:
                lons, lats = self.coordinates()
            egm96 = EGM96.EGM96reader()
            self.r += egm96(lons,lats)

        #Compatibility with Topo class
        lon, lat0 = self.coordinates(np.arange(self.ny),np.repeat((0,),self.ny))
        lon0, lat = self.coordinates(np.repeat((0,),self.nx),np.arange(self.nx))
        self.lon = lon
        self.lat = lat
        self.topo = self.r
        self.dlon = self.xres
        self.dlat = self.yres
        self.region = self.extent


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
            #singularities at edges
            distx[:,0] = distx[:,1]
            disty[0] = disty[1]
        else:
            distx = np.abs(self.xres)
            disty = np.abs(self.yres)

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



    def shaded_relief(self,azdeg=225,altdeg=30,cmap=pl.get_cmap('Greys'),alpha=1,aspect=None,smoothing=None,vmin='default',vmax='default',downsampl=1):
        """
        Function to plot a shaded relief of the DEM
        azdeg : f, azimuth (measured anti-clockwise from east) of the light source in degress for the shaded relief
        altdeg : f, altitude (measured up from the plane of the surface) of the light source in degrees
        cmap : object of class Colormap, default 'Greys'
        alpha : f, 0 to 1, transparency
        asp : f, imshow window aspect, default is 1
        smoothing : int, apply a gaussian filter of size 'smoothing'. Can be used to reduce noise (Default is none)
        vmin/vmax : f, setting these values manually can reduce areas of saturation in the plot, default is min/max of the DEM
        downsampl : int, factor to downsample the matrix if too large
        """

        # Smooth elevation
        if smoothing!=None:
            data = gaussian_filter(self.r,smoothing)
        else:
            data = self.r

        # down sample data
        if downsampl!=1:
            data = data[::downsampl,::downsampl]


        # mask potential Nan values
        data = np.ma.masked_array(data,mask=np.isnan(data))

        # set vmin/vmax value
        if vmin=='default':
            vmin = data.min()
        if vmax=='default':
            vmax=data.max()

        # compute shaded relief
        ls = LightSource(azdeg=azdeg, altdeg=altdeg)
        rgb0 = cmap((data - vmin) / (vmax - vmin))
        rgb1 = ls.shade_rgb(rgb0, elevation=data)

        # plot
        pl.imshow(rgb1,extent=self.extent,interpolation='bilinear',alpha=alpha,aspect=aspect)

        return rgb1

    
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
        
            

    def altitudinal_analysis(self,indata,alt_bins,operators=(np.mean,),srcnodata=np.nan,plot=False,ylim='default'):
        """
        Perform an altitude analysis of indata.
        For each altitude bin alt_bins, compute the statistics defined by the operators (default is np.mean)
        Optionally, display a plot of the distribution of mean and std as a function of altitude
        Arguments :
        - indata : np.array, data to be analysed, must have same shape as self.r
        - alt_bins : array, the altitude intervalls for the analysis
        - operators : list of operators to compute the statistics (e.g : np.mean, np.std, np.min...), default is np.mean
        - srcnodata : f, no data value for input
        - plot : bool, if set to True display a plot of the distribution of mean and std as a function of altitude
        - ylim : the y limits for the plot
        Outputs :
        - stats : np.array, (operators size) x (alt_bins size -1), contains each statistics for each altitude intervall
        """

        nz = len(alt_bins)-1
        nstats = len(operators)
        stats = np.nan*np.zeros((nstats,nz))

        for i in xrange(nz):
            data = indata[(self.r>alt_bins[i]) & (self.r<=alt_bins[i+1])]
            data = data[np.isfinite(data)];
            data = data[data!=srcnodata];

            #remove outliers
            std = np.std(data)
            mean = np.mean(data)
            data = data[(data>mean-3*std) & (data<mean+3*std)]

            #Compute statistics
            if len(data)!=0:
                for k in xrange(len(operators)):
                    stats[k,i] = operators[k](data)

        #Plot
        if plot==True:
            alt_mid = (alt_bins + np.roll(alt_bins,-1))/2  #Intervall midpoint
            x = alt_mid[:-1]
            xerr = x-alt_bins[:-1]   #Intervall half-width

            pl.errorbar(x,stats[0],xerr=xerr,yerr=stats[1],fmt='o',color='k') #display with error bars
            pl.hlines(0,np.min(x-xerr),np.max(x+xerr),linestyles='dashed')
            if ylim!='default':
                pl.ylim(ylim)

        return stats


    def along_slope_gradient(self,sigmax,sigmay):
        """
        Compute the DEM gradient (=slope) after smoothing the DEM along the slope with a gaussian kernel of size (sigmax, sigmay).
        Particularly useful to smooth the slope on glaciers without introducing to much errors from the edges.
        """

        # first smooth gradient, to determine slope direction
        dem_sm = gaussian_filter(self.r,1)
        xgrad, ygrad = np.gradient(dem_sm)
        alpha = np.arctan2(xgrad,ygrad)

        # smooth along slope
        zsm = np.nan*np.zeros_like(self.r)
        for i,j in np.transpose(np.where(np.isfinite(self.r))):
            f = gaussian_kernel(sigmax,sigmay,alpha[i,j])
            xsize, ysize = f.shape
            xmin = i-xsize/2
            xmax = i+xsize/2+1
            ymin = j-ysize/2
            ymax = j+ysize/2+1
            if ((xmin>0) & (xmax<self.r.shape[0]) & (ymin>0) & (ymax<self.r.shape[1])):
                data=self.r[xmin:xmax,ymin:ymax]
                nan=f[np.isfinite(data)]
                zsm[i,j] = np.nansum((data*f))/np.nansum(nan)
    
        
        # compute final slope
        if self.srs.IsProjected()==0:
            lon, lat = self.coordinates()
            lonr = np.roll(lon,1,1)
            latl = np.roll(lat,1,0)
            distx = geo.dist_ortho(lon,lat,lonr,lat)
            disty = geo.dist_ortho(lon,lat,lon,latl)
            #singularities at edges
            distx[:,0] = distx[:,1]
            disty[0] = disty[1]
        else:
            distx = np.abs(self.xres)
            disty = np.abs(self.yres)

        #Compute z gradients
        f1 = np.array([[-1,0,1],[-2,0,2],[-1,0,1]])
        f2 = f1.transpose()
        g1 = signal.convolve(zsm,f1,mode='same')
        g2 = signal.convolve(zsm,f2,mode='same')
        mpm = np.sqrt((g1/distx)**2 + (g2/disty)**2)/8
        slope = np.arctan(mpm)*180/np.pi
        
        return slope


def gaussian_kernel(sigmax,sigmay,alpha,fmin=0.01):
    """
    Return a normalized gaussian kernel with standard deviation "sigmax" and "sigmay" along the axis rotated by an angle "alpha" in radians. 
    The size of the matrix is calculated as 3*max(sigmax,sigmay).
    All elements below fmin are set to 0 and matrix is cut so that no line/row contains only zeros, thus matrix size is not always square.
    """

    size=np.ceil(max(sigmax,sigmay)*3)
    size=int(size//2*2+1)  # make it an odd number
    #if size%2!=1:
    #    print "size must be an odd number"
    #    return 0

    # create meshgrid
    x = np.arange(-size/2,size/2)+1
    y = np.arange(-size/2,size/2)+1
    xx, yy = np.meshgrid(x,y)

    # rotate axes
    xp = xx*np.cos(alpha)+yy*np.sin(alpha)
    yp = -xx*np.sin(alpha)+yy*np.cos(alpha)

    #compute gaussian kernel
    f = np.exp(-(xp**2/(2*sigmax**2)+yp**2/(2*sigmay**2)))
    f = f/np.sum(f)

    # cut the matrix to values higher than 0.01
    f=f[:,~np.all(f<0.01, axis=0)]
    f=f[~np.all(f<0.01, axis=1),:]
    f = f/np.sum(f)

    return f



"""
Additional shaded relief algorithm
Copied from http://rnovitsky.blogspot.co.uk/2010/04/using-hillshade-image-as-intensity.html
Example : 
> rgb=set_shade(input,cmap=plt.get_cmap('gist_earth'),scale=1.0)
> plt.imshow(rgb)
> plt.show()


def set_shade(a,intensity=None,cmap=pl.get_cmap('jet'),scale=10.0,azdeg=165.0,altdeg=45.0):
    ''' sets shading for data array based on intensity layer
    or the data's value itself.
    inputs:
    a - a 2-d array or masked array
    intensity - a 2-d array of same size as a (no chack on that)
    representing the intensity layer. if none is given
    the data itself is used after getting the hillshade values
    see hillshade for more details.
    cmap - a colormap (e.g matplotlib.colors.LinearSegmentedColormap
    instance)
    scale,azdeg,altdeg - parameters for hilshade function see there for
    more details
    output:
    rgb - an rgb set of the Pegtop soft light composition of the data and 
    intensity can be used as input for imshow()
    based on ImageMagick's Pegtop_light:
    http://www.imagemagick.org/Usage/compose/#pegtoplight'''
    if intensity is None:
        # hilshading the data
        intensity = hillshade(a,scale=scale,azdeg=165.0,altdeg=45.0)
    else:
        # or normalize the intensity
        intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
        # get rgb of normalized data based on cmap

    rgb = cmap((a-a.min())/float(a.max()-a.min()))[:,:,:3]
    # form an rgb eqvivalent of intensity
    d = intensity.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2*d*rgb+(rgb**2)*(1-2*d)
    
    return rgb


def hillshade(data,scale=10.0,azdeg=165.0,altdeg=45.0):
    ''' convert data to hillshade based on matplotlib.colors.LightSource class.
    input:
    data - a 2-d array of data
    scale - scaling value of the data. higher number = lower gradient
    azdeg - where the light comes from: 0 south ; 90 east ; 180 north ;
    270 west
    altdeg - where the light comes from: 0 horison ; 90 zenith
    output: a 2-d array of normalized hilshade
    '''
    # convert alt, az to radians
    az = azdeg*np.pi/180.0
    alt = altdeg*np.pi/180.0
    # gradient in x and y directions
    dx, dy = np.gradient(data/float(scale))
    slope = 0.5*np.pi - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dx, dy)
    intensity = np.sin(alt)*np.sin(slope) + np.cos(alt)*np.cos(slope)*np.cos(-az - aspect - 0.5*np.pi)
    intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
    return intensity
"""
