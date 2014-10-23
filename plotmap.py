#coding=utf-8
"""
plotmap.py

Map : class to generate Basemap figures.

Creates a geo-referenced Basemap figure, and provides a number of methods 
to add different types of layers and information.

The methods are broadly arranged in the expected calling sequence. See 
method docstrings for further information.

For a practical example look at landsat_velocities/PlotFunctions.plot_map().

Author : Andrew Tedstone
Date: October 2014
"""

import numpy as np
import pylab as pl
from matplotlib import cm
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
from matplotlib.colors import LightSource

import georaster


class Map:

    map = None
    fig = None
    extent = None

    def __init__(self,ds_file=None,extent=None,lon_0=None):
        """

        Create the Map object, which itself creates a Basemap instance.

        The Map object must be initialised with georeferencing information.
        This can be provided in two ways:
            ds_file : str, link to a GeoTIFF. The extent of the plotting area
                      and the lon_0 will be set according to the properties of
                      the GeoTIFF.
        Or: 
            extent : (lon_lower_left,lon_upper_right,lat_lower_left,
                      lat_upper_right)
            lon_0 : float, longitude of origin

        E.g.
        >>> mymap = Map(ds_file='myim.tif')

        E.g.
        >>> mymap = Map(extent=(-50,-48,67,68),lon_0=70)

        """

        # Create figure
        self.fig = pl.figure()

        # Get basic georeferencing info for map
        # From a geoTIFF
        if ds_file <> None:
            ds = georaster.SingleBandRaster(ds_file,load_data=False)
            extent = ds.get_extent_latlon()  
            lon_0 = ds.srs.GetProjParm('central_meridian')
        # Otherwise check that it has been provided manually
        else:
            if (extent == None) or (lon_0 == None):
                print 'Either ds_file must be provided, or extent and lon_0.'
                raise AttributeError

        self.extent = extent
        lonll, lonur, latll, latur = extent

        # Create Basemap    
        self.map = Basemap(llcrnrlon=lonll,llcrnrlat=latll,urcrnrlon=lonur,
                      urcrnrlat=latur,resolution='i',projection='tmerc',
                      lon_0=lon_0,lat_0=0)



    def plot_background(self,bg_file,region='all',coarse=False):
        """
        Plot a background image onto the Basemap, in black and white.

        Optionally, coarsen resolution of background image, which will
        result in smaller file sizes in saved vector formats.

        Inputs:
            bg_file : str, path and filename of single-band GeoTIFF
            region : 'all' or latlon tuple (lonll,lonur,latll,latur)
            coarse : False or int to coarsen image by (e.g. 2)

        """

        bg_im = georaster.SingleBandRaster(bg_file,load_data=False)
        if region <> 'all':
            try:
                bg_im.r = bg_im.read_single_band_subset(region,latlon=True) 
            except ValueError:
                print 'Region not within background image bounds, using full image bounds instead . . .'
                bg_im.r = bg_im.read_single_band()    
        else:
            bg_im.r = bg_im.read_single_band()

        # Reduce image resolution
        if coarse <> False:
            if type(coarse) == int:
                bg_im.r = bg_im.r[::coarse,::coarse]

        bg_im.r = np.where(bg_im.r == 0,np.nan,bg_im.r)   #remove black color
        pl.imshow(bg_im.r,cmap=cm.Greys_r,
                  extent=bg_im.get_extent_projected(self.map),
                  interpolation='nearest')



    def plot_dem(self,dem_file,region='all'):
        """
        Plot a DEM using light-shading on the Basemap.

        Inputs:
            dem_file : path and filename of GeoTIFF DEM.
            region : 'all' or latlon tuple (lonll,lonur,latll,latur)

        """

        dem_im = georaster.SingleBandRaster(dem_file,load_data=False)   
        if region <> 'all':
            try:
                dem_im.r = dem_im.read_single_band_subset(region,latlon=True) 
            except ValueError:
                print 'Region not within DEM bounds, using full DEM bounds instead . . .'
                dem_im.r = dem_im.read_single_band()   
        else:
            dem_im.r = dem_im.read_single_band()
        ls = LightSource(azdeg=135,altdeg=80)
        rgb = ls.shade(dem_im.r,cmap=cm.Greys_r)  
        pl.imshow(rgb,extent=dem_im.get_extent_projected(self.map),
            interpolation='nearest')



    def plot_mask(self,mask_file,color='turquoise',region='all'):
        """
        Plot masked values of the provided dataset.

        This can be useful to display 'bad' regions. If the colormap is set to
        discrete, the regions will be plotted in red, otherwise in white.

        Arguments:
            mask_file : str, path to geoTIFF of mask
            data : georaster object of the data to mask
            color : any matplotlib color, str or RGB triplet
            region : optional, latlon tuple (lonll,lonur,latll,latur) 

        """
            
        mask = georaster.SingleBandRaster(mask_file,load_data=False)
        if region == 'all':
            mask.r = mask.read_single_band_subset(region,latlon=True)
        else:
            mask.r = mask.read_single_band()

        #Pixels outside mask are transparent
        mask.r = np.where(mask.r==0,np.nan,1)

        # Now plot the bad data values
        cmap = cm.jet
        if color == 'turquoise':
            cmap.set_over((64./255,224./255,208./255))
            alpha=0.6
        else:
            try:
                cmap.set_over(eval(color))  #color is a RGB triplet
            except NameError: 
                cmap.set_over(color)      #color is str, e.g red, white...
            alpha=1

        pl.imshow(mask.r,extent=mask.get_extent_projected(self.map),
                  cmap=cmap,vmin=-1,vmax=0,interpolation='nearest',alpha=alpha)



    def plot_data(self,data,vmin='min',vmax='max',cmap='cm.jet'):
        """
        Basic function to plot a dataset with minimum and maximum values.

        In many cases you will be better off calling pl.imshow yourself
        instead.

        Arguments:
            data : georaster object of data to plot
            vmin : 'min' or minimum value to plot
            vmax : 'max' or maximum value to plot
            cmap : colormap, in format cn.<name>

        """

        if vmin == 'min':
            vmin = np.nanmin(data.r)
        if vmax == 'max':
            vmax = np.nanmax(data.r)
        # Discrete colormap option
        if cmap=='discr':
            cmap = cm.Blues
            bounds = [0,5,10,15,20,30,40,80,120,200]
            norm = colors.BoundaryNorm(bounds, cmap.N)
            pl.imshow(data.r,extent=data.get_extent_projected(self.map),
                cmap=cmap,
                vmin=0,norm=norm,interpolation='nearest',alpha=1)
        # Continuous colormap option
        else:
            cmap = eval(cmap)
            pl.imshow(data.r,extent=data.get_extent_projected(self.map),
                cmap=cmap,
                vmin=vmin,vmax=vmax,interpolation='nearest',alpha=1)


    def geo_ticks(self,mstep,pstep):
        """
        Add geographic (latlon) ticks to plot.

        Arguments:
            mstep :  float, meridians stepping, degrees
            pstep : float, parallels stepping, degrees

        """
        lonll, lonur, latll, latur = self.extent
        m0 = int(lonll/mstep)*mstep
        m1 = int(lonur/mstep+1)*mstep
        p0 = int(latll/pstep)*pstep
        p1 = int(latur/pstep+1)*pstep
        self.map.drawparallels(np.arange(p0,p1,pstep),labels=[1,0,0,0],
            linewidth=0)
        self.map.drawmeridians(np.arange(m0,m1,mstep),labels=[0,0,0,1],
            linewidth=0)



    def plot_scale(self,length):
        """
        Plot a scale bar on the figure.

        Arguments:
            length : float, length of scale bar in map units/

        """
        lonll, lonur, latll, latur = self.extent
        xloc = lonll + 0.85*(lonur-lonll)
        yloc = latll + 0.93*(latur-latll)
        self.map.drawmapscale(xloc,yloc,(lonur+lonll)/2,(latur+latll)/2,length,
                         barstyle='fancy')



    def plot_colorbar(self,label=None):
        """
        Draw colorbar of the same height as the figure.

        Arguments:
            label : str, text to label colorbar with

        """
        ax = self.fig.gca()   #get axis properties
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = pl.colorbar(cax=cax)
        if label != None:
            cb.set_label(label)
        cb.set_alpha(1)  #no transparency
        cb.draw_all()



    def save_figure(self,outfile):
        """
        Save figure to outfile, having adjusted subplots, with 300dpi resn.

        Arguments:
            outfile : str, path and filename of file to save to.

        """

        self.fig.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.07)
        self.fig.savefig(outfile,dpi=300)
