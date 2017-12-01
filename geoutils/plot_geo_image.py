#!/usr/bin/env python
#coding=utf-8

"""
Description : Visualisation of any image that is compatible with GDAL

Author : Amaury Dehecq
Date : 11/2017
"""


#Python libraries
from __future__ import print_function
import os, sys, argparse
import matplotlib.pyplot as plt
import numpy as np

#Personal libraries
import georaster as raster




if __name__=='__main__':

    #Set up arguments
    parser = argparse.ArgumentParser(description='Visualisation tool for any image supported by GDAL.')

    #Positional arguments
    parser.add_argument('filename', type=str, help='str, path to the image')

    #optional arguments
    parser.add_argument('-cmap', dest='cmap', type=str, default='default', help='str, a matplotlib colormap string (default is from rcParams).')
    parser.add_argument('-vmin', dest='vmin', type=str, default='default', help='float, the minimum value for colorscale (default is calculated min value).')
    parser.add_argument('-vmax', dest='vmax', type=str, default='default', help='float, the maximum value for colorscale (default is calculated max value).')
    parser.add_argument('-band', dest='band', type=int, default=0, help='int, which band to display (start at 0) for multiband images (Default is 0).')
    parser.add_argument('-nocb', dest='nocb', help='If set, will not display a colorbar (Default is to display the colorbar).',action='store_true')
    parser.add_argument('-clabel', dest='clabel', type=str, default='', help='str, the label for the colorscale (Default is empty).')
    parser.add_argument('-title', dest='title', type=str, default='', help='str, figure title (Default is empty).')
    parser.add_argument('-figsize', dest='figsize', type=str, default='default', help='str, figure size, must be a tuple of size 2, either written with quotes, or two numbers seperated by coma, no space (Default is from rcParams).')
    parser.add_argument('-max_size', dest='max_size', type=int, default=2000, help='int, image size is limited to max_size**2 for faster reading/displaying (Default is 2000).')
    parser.add_argument('-save', dest='save', type=str, default='', help='str, filename to the output filename to save to disk (Default is displayed on screen).')
    parser.add_argument('-dpi', dest='dpi', type=str, default='default', help='int, dpi value to use when saving figure (Default is from rcParams).')
    parser.add_argument('-nodata', dest='nodata', type=str, default='default', help='float, no data value (Default is read from file metadata).')
    
    args = parser.parse_args()

    ## Load image metadata ##
    img = raster.MultiBandRaster(args.filename,load_data=False)
    xmin, xmax, ymin, ymax = img.extent

    if args.nodata=='default':
        nodata = img.ds.GetRasterBand(1).GetNoDataValue()
    else:
        try:
            nodata = float(args.nodata)
        except ValueError:
            print("ERROR: nodata must be a float, currently set to %s" %args.nodata)
            sys.exit(1)


    ## Resample if image is too large ##
    if img.nx*img.ny>args.max_size**2:
        step = max(int(img.nx/args.max_size),int(img.ny/args.max_size))
        print("Image will be downsampled by a factor %i." %step)
    else:
        step=1

    ## Read image ##
    img = raster.MultiBandRaster(args.filename,downsampl=step,bands=(args.band+1,))
    if nodata!=None:
        data = np.ma.masked_array(img.r[:,:,0],mask=(img.r[:,:,0]==nodata))
    else:
        data = img.r[:,:,0]

    ## Set default parameters ##

    # vmin
    if args.vmin=='default':
        vmin=np.nanmin(data)
    else:
        try:
            vmin=float(args.vmin)
        except ValueError:
            print("ERROR: vmin must be a float, currently set to %s" %args.vmin)
            sys.exit(1)

    # vmax
    if args.vmax=='default':
        vmax=np.nanmax(data)
    else:
        try:
            vmax=float(args.vmax)
        except ValueError:
            print("ERROR: vmax must be a float, currently set to %s" %args.vmax)
            sys.exit(1)

    # color map
    if args.cmap=='default':
        cmap = plt.rcParams['image.cmap']
    elif args.cmap in plt.cm.datad.keys():
        cmap = args.cmap
    else:
        print("ERROR: cmap set to %s, must be in:" %args.cmap)
        for i, elem in enumerate(plt.cm.datad.keys(), 1):
            print(str(elem), end='\n' if i % 10 == 0 else ', ')
        sys.exit(1)

    # Figsize
    if args.figsize=='default':
        figsize = plt.rcParams['figure.figsize']
    else:
        print(eval(args.figsize))
        print(tuple(eval(args.figsize)))
        try:
            figsize=tuple(eval(args.figsize))
            xfigsize, yfigsize = figsize
        except:
            print("ERROR: figsize must be a tuple of size 2, currently set to %s" %args.figsize)
            sys.exit(1)

    # dpi
    if args.dpi=='default':
        dpi = plt.rcParams['figure.dpi']
    else:
        try:
            dpi = int(args.dpi)
        except ValueError:
            print("ERROR: dpi must be an integer, currently set to %s" %args.dpi)
            sys.exit(1)


    
    ## Plot data ##
    
    fig = plt.figure(figsize=figsize)

    # plot
    plt.imshow(data,extent=(xmin,xmax,ymin,ymax),cmap=cmap,interpolation='nearest',vmin=vmin,vmax=vmax)

    # colorbar
    if args.nocb==False:
        cb = plt.colorbar()

        if args.clabel!='':
            cb.set_label(args.clabel)

    # title
    if args.title!='':
        plt.title(args.title)

    plt.tight_layout()

    # Save
    if args.save!='':
        plt.savefig(args.save,dpi=dpi)
        print("Figure saved to file %s." %args.save)
    else:
        print("Figure displayed on screen.")
        plt.show()
        
