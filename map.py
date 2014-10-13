#coding=utf-8
"""
Class to generate Basemap figures.

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

        Provide either: 
            ds_file (str, link to geoTIFF)
        Or: 
            extent : (lonll,lonur,latll,latur)
            lon_0 : float, longitude of origin

        """

        # Create figure
        self.fig = pl.figure()
        self.plt = pl.subplot(111)

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


    def plot_background(self,bg_file,region='all'):
        bg_im = georaster.SingleBandRaster(bg_file,read_data=False)
        if region <> 'all':
            bg_im.r = bg_im.read_single_band_subset(region,latlon=True) 
        else:
            bg_im.r = bg_im.read_single_band()
        # Reduce image resolution
        bg_im.r = bg_im.r[::4,::4]
        bg_im.r = np.where(bg_im.r == 0,np.nan,bg_im)   #remove black color
        self.plt.imshow(bg_im.r,cmap=cm.Greys_r,extent=bg_im.get_extent_projected(map),
                  interpolation='nearest')


    def plot_dem(self,dem_file,region='all'):

        dem_im = georaster.SingleBandRaster(dem_file,load_data=False)   
        if region <> 'all':
            dem_im.r = dem_im.read_single_band_subset(region,latlon=True) 
        else:
            dem_im.r = dem_im.read_single_band()
        ls = LightSource(azdeg=135,altdeg=80)
        rgb = ls.shade(dem_im.r,cmap=cm.Greys_r)  
        self.plt.imshow(rgb,extent=dem_im.get_extent_projected(self.map),interpolation='nearest')




    def plot_mask(self,mask_file,data,region='all'):
        """

        Arguments:
            mask_file : str, path to geoTIFF of mask
            data : georaster object of the data to mask
            region : optional, (lonll,lonur,latll,latur) of area to plot

        """
            
        mask = georaster.SingleBandRaster(maskfile,load_data=False)
        if region == 'all':
            mask.r = mask.read_single_band_subset(region,latlon=True)
        else:
            mask.r = mask.read_single_band()

        # set NaN value to -999 to display in a different color
        data.r = np.where(np.isnan(data.r),-999,data.r)  
        
        data_nan = np.where((mask.r > 0) & (data.r == -999),data.r,np.nan)
        data.r = np.where((mask.r > 0) & (data.r > 0),data.r,np.nan)

        # Now plot the bad data values
        red_cmap = cm.jet
        if cmap == 'discr':
            red_cmap.set_under((155./255,0./255,0./255)) #gaps displayed in red
        else:
            red_cmap.set_under((1,1,1))  #gaps displayed in white
        self.plt.imshow(data_nan,extent=data.get_extent_projected(self.map),
                  cmap=red_cmap,vmin=-1,vmax=0,interpolation='nearest',alpha=1)


    def plot_data(self,data,vmin='min',vmax='max',cmap='cm.jet'):
        """

        Arguments:
            data : georaster object of data to plot

        """
        print data.r

        if vmin == 'min':
            vmin = np.nanmin(data.r)
        if vmax == 'max':
            vmax = np.nanmax(data.r)
        # Discrete colormap option
        if cmap=='discr':
            cmap = cm.Blues
            bounds = [0,5,10,15,20,30,40,80,120,200]
            norm = colors.BoundaryNorm(bounds, cmap.N)
            self.plt.imshow(data.r,extent=data.get_extent_projected(self.map),cmap=cmap,
                vmin=0,norm=norm,interpolation='nearest',alpha=alpha)
        # Continuous colormap option
        else:
            cmap = eval(cmap)
            self.plt.imshow(data.r,extent=data.get_extent_projected(self.map),cmap=cmap,
                vmin=vmin,vmax=vmax,interpolation='nearest',alpha=1)


    def plot_scales(self,mstep,pstep,length):
        lonll, lonur, latll, latur = self.extent
        m0 = int(lonll/mstep)*mstep
        m1 = int(lonur/mstep+1)*mstep
        p0 = int(latll/pstep)*pstep
        p1 = int(latur/pstep+1)*pstep
        xloc = lonll + 0.85*(lonur-lonll)
        yloc = latll + 0.93*(latur-latll)
        self.map.drawparallels(np.arange(p0,p1,pstep),labels=[1,0,0,0],linewidth=0)
        self.map.drawmeridians(np.arange(m0,m1,mstep),labels=[0,0,0,1],linewidth=0)
        self.map.drawmapscale(xloc,yloc,(lonur+lonll)/2,(latur+latll)/2,length,
                         barstyle='fancy')


    def plot_colorbar(self):
        # Draw colorbar of the same height as the figure
        ax = self.plt   #get axis properties
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = pl.colorbar(cax=cax)
        if label!='':
            cb.set_label(label)
        cb.set_alpha(1)  #no transparency
        cb.draw_all()


    def save_figure(self,outfile):
        self.fig.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.07)
        self.fig.savefig(outfile,dpi=300)
