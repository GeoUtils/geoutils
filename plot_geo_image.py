#!/usr/bin/python
#coding=utf-8

"""
Description : Visualisation of any image that is compatible with GDAL

Author : Amaury Dehecq
Date : 02/2014
"""


#Python libraries
import os, sys
import pylab as pl
from matplotlib import cm
import numpy as np

#Personal libraries
import georaster as geo




if __name__=='__main__':

    # =============================================================================
    def Usage():
        print '*** Plot any image supported by GDAL ***'
        print
        print 'Usage : plot_raster.py filename colorbar=False cmap=jet'
        print
        print 'input parameters :'
        print '  filename : str, path to the image file'
        print 
        print 'Opional arguments :'
        print '  colorbar : True or False, display colorbar (Default is False)'
        print '  cmap : matplotlib colormap (default is jet)'
        print 
# =============================================================================
    
    if len(sys.argv)<2:
        Usage()
        sys.exit(1)

    filename = sys.argv[1]
    
    #read optional arguments
    kwargs = {}
    items = ['colorbar','cmap','vmin','vmax']

    kwargs['colorbar']=False
    kwargs['cmap']='jet'

    if len(sys.argv)>2:
        for arg in sys.argv[2:]:
    
            item, value = arg.split('=')
    
            if item in items:
                kwargs[item]=value

    for key in kwargs.keys():
        if kwargs[key]=='True':
            kwargs[key] = True
        elif kwargs[key]=='False':
            kwargs[key] = False


    img = geo.SingleBandRaster(filename)
    data = img.r
    xmin, xmax, ymin, ymax = img.extent

    cmap = eval('cm.'+kwargs['cmap'])

    #Resample if image is too large
    ysize, xsize = data.shape
    if xsize*ysize>2000*2000:
        step = max(int(xsize/2000),int(ysize/2000))
        data = data[::step,::step]

    if kwargs.has_key('vmin')==False:
        kwargs['vmin']=np.nanmin(data)
    else:
        kwargs['vmin']=float(kwargs['vmin'])

    if kwargs.has_key('vmax')==False:
        kwargs['vmax']=np.nanmax(data)
    else:
        kwargs['vmax']=float(kwargs['vmax'])

    pl.imshow(data,extent=(xmin,xmax,ymin,ymax),cmap=cmap,interpolation='none',vmin=kwargs['vmin'],vmax=kwargs['vmax'])
    if kwargs['colorbar']==True:
        pl.colorbar()
    pl.show()
