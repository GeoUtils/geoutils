#!/usr/bin/env python
#coding=utf-8

"""
Description : Calculate the difference between two georeferenced rasters.

Author : Amaury Dehecq
"""

#Python libraries
import argparse, os
import numpy as np
import gdal

#Personal libraries
import georaster as geor

## Dictionary to convert between resampling algo names and GDAL codes ##
GRA_codes = {'near':gdal.GRA_NearestNeighbour, 'bilinear':gdal.GRA_Bilinear, 'cubic':gdal.GRA_Cubic, 'cubicspline':gdal.GRA_CubicSpline, 'lanczos':gdal.GRA_Lanczos, 'average':gdal.GRA_Average, 'mode':gdal.GRA_Mode,  'max':gdal.GRA_Max, 'min':gdal.GRA_Min, 'med':gdal.GRA_Med, 'Q1':gdal.GRA_Q1, 'Q3':gdal.GRA_Q3}



def geodiff(im1,im2,imout,inverse=False,resampling=1):
    """
    Calculate the difference between images im1 and im2. First, im1 is loaded in area of overlap with im2 and im2 is reprojected on im1's grid. Then the difference im1-im2 is calculated. If inverse is set to True, im2-im1 is calculated instead. 
    Inputs:
    im1: str, path to the first image
    im2: str, path to the second image
    imout: str, path to the output image
    inverse: optional - bool, set to True to calculate im2-im1 (Default is False)
    resampling: optional - int, resampling algorithm. Use the gdal.GRA codes: 0 is nearest, 1 is bilinear etc. (Default is 1 - bilinear)
    """

    # Read first image
    img1 = geor.SingleBandRaster(im1, load_data=False)
    xmin, xmax, ymin, ymax = img1.intersection(im2)
    #xmin, xmax, ymin, ymax = img1.extent
    nodata1 = img1.ds.GetRasterBand(1).GetNoDataValue()
    
    # Load second image
    img2 = geor.SingleBandRaster(im2,load_data=False)
    dtype2 = img2.ds.GetRasterBand(1).DataType
    nodata2 = img2.ds.GetRasterBand(1).GetNoDataValue()

    # Read first image
    img1 = geor.SingleBandRaster(im1, load_data=(xmin,xmax,ymin,ymax),latlon=False)
    
    # Reproject img2 if needed
    if ((img1.nx != img2.nx) or (img1.ny != img2.ny)):
        img2_proj = img2.reproject(img1.srs, nx=img1.nx, ny=img1.ny, xmin=xmin, ymax=ymax, xres=img1.xres, yres=img1.yres, dtype=dtype2, nodata=nodata2, interp_type=resampling, progress=True)
    else:
        img2_proj = geor.SingleBandRaster(im2)

    # Calculate difference
    if inverse==False:
        diff = img1.r - img2_proj.r
    elif inverse==True:
        diff = img2_proj.r - img1.r
    else:
        raise ValueError("'inverse' must be True or False")

    out_nodata = -32767 #np.finfo('float').min
    diff[img1.r==nodata1] = out_nodata
    diff[img2_proj.r==nodata2] = out_nodata
    
    # Save to GTiff
    return geor.simple_write_geotiff(imout, diff, img1.ds.GetGeoTransform(), wkt=img1.srs.ExportToWkt(), dtype=6, nodata_value=out_nodata, options=None)

    
if __name__=='__main__':

    #Set up arguments
    parser = argparse.ArgumentParser(description="Calculate the difference between images im1 and im2.  First, im1 is loaded in area of overlap with im2 and im2 is reprojected on im1's grid. Then the difference im1-im2 is calculated. If inverse is set to True, im2-im1 is calculated instead.")
    
    #Positional arguments
    parser.add_argument('im1', type=str, help='str, path to the first image')
    parser.add_argument('im2', type=str, help='str, path to the second image')
    parser.add_argument('imout', type=str, help='str, path to the output file')

    # Optional arguments
    parser.add_argument('-i', dest='inverse', action='store_true', help='If set, will calculate im2-im1.')
    parser.add_argument('-r', dest='resampling', type=str, default='bilinear', help="str, GDAL resampling algorithms used to reproject im2 onto im1's grid (Default is 'near')")
    args = parser.parse_args()


    geodiff(args.im1,args.im2,args.imout,inverse=args.inverse, resampling=GRA_codes[args.resampling])
    
