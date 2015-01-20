#!/usr/bin/python
#coding=utf-8

"""
Description : Crop and project an image in the same grid as a reference image

Author : Amaury Dehecq
"""

#Python libraries
from osgeo import gdal, osr
import sys, os

#Personal libraries
import georaster as geo


# =============================================================================
def Usage():
    print '*** Crop and project an image in the same grid as a reference image ***'
    print
    print 'Usage : crop2image.py im_ref im2crop outfile'
    print
    print 'input parameters :'
    print '  im_ref : str, path to the reference image'
    print '  im2crop : str, path to the image to crop'
    print '  outfile : str, name of the output file'
    print 
# =============================================================================
    

if len(sys.argv)<4:
    Usage()
    sys.exit(1)


#Read arguments
im_ref = sys.argv[1]
im2crop = sys.argv[2]
outfile = sys.argv[3]


#Read reference image spatial reference system
ds = gdal.Open(im_ref)
srs_dest=osr.SpatialReference()
wkt = ds.GetProjectionRef()
srs_dest.ImportFromWkt(wkt)
proj=srs_dest.ExportToProj4()


#Read reference image size and corners coordinates
img = geo.SingleBandRaster(im_ref)
pixelWidth, pixelHeight = img.get_pixel_size()
xmin, xmax, ymin, ymax = img.extent


#Crop and reproject
cmd = "gdalwarp -te %.8f %.8f %.8f %.8f -tr %.8f %.8f -t_srs '%s' %s %s -overwrite" %(xmin,ymin,xmax,ymax,pixelWidth,pixelHeight,proj,im2crop,outfile)
print cmd; os.system(cmd)
