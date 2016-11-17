#!/usr/bin/env python
#coding=utf-8

"""
Description : Crop two images to the common extent. Convenient for feature-tracking.

Author : Amaury Dehecq
"""

#Python libraries
import argparse, os

#Personal libraries
import georaster as raster



def crop2common(im1,im2,im1out,im2out):
    """
    Crop two images to the common extent.
    Inputs:
    im1: str, path to the first image
    im2: str, path to the second image
    im1out: str, path to the cropped first image
    im2out: str, path to the cropped second image
    """

    #Read first image extent
    img1 = raster.SingleBandRaster(im1,load_data=False)
    xmin1, xmax1, ymin1, ymax1 = img1.extent

    #Read first image extent
    img2 = raster.SingleBandRaster(im2,load_data=False)
    xmin2, xmax2, ymin2, ymax2 = img2.extent

    # compute common extent
    xmin = max(xmin1,xmin2)
    xmax = min(xmax1,xmax2)
    ymin = max(ymin1,ymin2)
    ymax = min(ymax1,ymax2)
                
    #Crop first image
    cmd = "gdalwarp -te %.8f %.8f %.8f %.8f %s %s -overwrite" %(xmin,ymin,xmax,ymax,im1,im1out)
    print cmd; os.system(cmd)

    #Crop second image
    cmd = "gdalwarp -te %.8f %.8f %.8f %.8f %s %s -overwrite" %(xmin,ymin,xmax,ymax,im2,im2out)
    print cmd; os.system(cmd)

    
if __name__=='__main__':

    #Set up arguments
    parser = argparse.ArgumentParser(description='Crop two images to the common extent')
    
    #Positional arguments
    parser.add_argument('im1', type=str, help='str,path to the first image')
    parser.add_argument('im2', type=str, help='str,path to the second image')
    parser.add_argument('im1out', type=str, help='str,path to the cropped first image')
    parser.add_argument('im2out', type=str, help='str,path to the second cropped image')
    args = parser.parse_args()


    crop2common(args.im1,args.im2,args.im1out,args.im2out)
    
