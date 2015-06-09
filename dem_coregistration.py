#!/usr/bin/env python
#coding=utf-8

"""
Description : Fine coregistration of 2 DEMs using the method presented in Nuth & Kaab 2011

Author : Amaury Dehecq
Date : June 2015
"""

#Python libraries
from scipy import ndimage
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as pl
from scipy.interpolate import RectBivariateSpline
import argparse

#Personal libraries
import georaster as raster

#Disable warnings
import warnings
warnings.filterwarnings("ignore")


def grad2d(dem):
  '''
  Calculate the slope and gradient of a DEM
  '''

#  I = ndimage.gaussian_filter(dem,0.333)
  g1, g2 = np.gradient(dem)

  slope_pix = np.sqrt(g1**2 + g2**2)
  aspect = np.arctan2(g2,g1)
  aspect = aspect+np.pi

  return slope_pix,aspect


def coregistration(master_dem,slave_dem,plot=False):
    """
    Compute the horizontal shift between 2 DEMs using the method presented in Nuth & Kaab 2011
    Inputs :
    - master_dem and slave_dem are the elevation matrices associated to each DEMs, matrices must have the same size and master_dem is used as the reference for geocoding
    """

    #Elevation difference
    dh = master_dem - slave_dem

    #compute terrain aspect
    slope_pix, aspect = grad2d(master_dem)
    slope_pix = np.where(np.isnan(dh),np.nan,slope_pix)
    aspect = np.where(np.isnan(dh),np.nan,aspect)
    target = dh/slope_pix


    #mean by slice
    slice_bounds = np.arange(0,2*np.pi,np.pi/36)
    mean=np.zeros([len(slice_bounds)])
    x_s=np.zeros([len(slice_bounds)])
    j=0
    for i in slice_bounds:
        target_slice = target[(i<aspect) & (aspect<i+np.pi/36)] #select target in the slice 
        target_slice = target_slice[(target_slice<200) & (target_slice>-200)] #avoid target>200 and target<-200
        mean[j] = np.mean(target_slice) #derive mean of target in the slice
        x_s[j] = i
        j=j+1

    #function to fit according to Nuth & Kaab
    x=aspect.ravel()
    y_meas = target.ravel()

    #remove non-finite values
    xf = x[(np.isfinite(x)) & (np.isfinite(y_meas))]
    yf = y_meas[(np.isfinite(x)) & (np.isfinite(y_meas))]

    #remove outliers
    p1 = np.percentile(yf,1)
    p99 = np.percentile(yf,99)
    xf = xf[(p1<yf) & (yf<p99)]
    yf = yf[(p1<yf) & (yf<p99)]


    #First guess
    p0 = (3*np.std(yf)/(2**0.5),0,np.mean(yf))

    #Least square fit   
    def peval(x,p):
        return p[0]*np.cos(p[1]-x) + p[2]

    def residuals(p,y,x):
        err = peval(x,p)-y
        return err

    plsq = leastsq(residuals, p0, args = (mean,x_s),full_output = 1)
    yfit = peval(x_s,plsq[0])

    #plotting results
    if plot==True:
        pl.plot(x_s,mean,'b.')
        pl.plot(x_s,yfit,'k-')
        #ax.set_ylim([np.min(mean),])
        pl.xlabel('Terrain aspect (rad)')
        pl.ylabel(r'dh/tan($\alpha$)')
        pl.show()

    a,b,c = plsq[0]
    east = a*np.cos(b)
    north = a*np.sin(b)

    return east, north, c

if __name__=='__main__':

    #Set up arguments
    parser = argparse.ArgumentParser(description='Fine coregistration of 2 DEMs using the method presented in Nuth & Kaab 2011')

    #Positional arguments
    parser.add_argument('master_dem', type=str, help='str,path to the master DEM')
    parser.add_argument('slave_dem', type=str, help='str, path to the slave DEM')
    parser.add_argument('outfile', type=str, help='str, path to the output coregistered DEM')

    #optional arguments
    parser.add_argument('-iter', dest='niter', type=int, default=5, help='int, number of iterations (default: 5)')
    parser.add_argument('-plot', dest='plot', help='Plot processing steps and final results',action='store_true')
    parser.add_argument('-m', dest='maskfile', type=str, default='none', help='str, path to a mask of same size as the master DEM, to filter out non stable areas such as glaciers (default is none)')
    parser.add_argument('-n1', dest='nodata1', type=str, default='none', help='int, no data value for master DEM if not specified in the raster file (default read in the raster file)')


    args = parser.parse_args()


    #Read DEMs
    master_dem = raster.SingleBandRaster(args.master_dem)
    master_dem.r = np.float32(master_dem.r)
    if args.nodata1!='none':
        master_dem.r[master_dem.r==int(args.nodata1)] = np.nan
    else:
        band=master_dem.ds.GetRasterBand(1)
        nodata = band.GetNoDataValue()
        master_dem.r[master_dem.r==nodata] = np.nan

    slave_dem = raster.SingleBandRaster(args.slave_dem)
    band=slave_dem.ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    slave_dem.r = np.float32(slave_dem.r)
    slave_dem.r[slave_dem.r==nodata] = np.nan

    #reproject slave DEM into the master DEM spatial reference system
    if master_dem.r.shape!=slave_dem.r.shape:
        band=master_dem.ds.GetRasterBand(1)
        dem2coreg = slave_dem.reproject(master_dem.srs, master_dem.nx, master_dem.ny, master_dem.extent[0], master_dem.extent[3], master_dem.xres, master_dem.yres, dtype=band.DataType, nodata=nodata, interp_type=1)

    else:
        dem2coreg = slave_dem.r

    dem2coreg[dem2coreg==nodata] = np.nan

    #fill NaN values for interpolation
    nanval = np.isnan(dem2coreg)
    slave_filled = np.where(np.isnan(dem2coreg),-9999,dem2coreg)
    
    #mask points
    if args.maskfile!='none':
        mask = raster.SingleBandRaster(args.maskfile)
        master_dem.r[mask.r>1] = np.nan


    #Set master DEM grid for later resampling
    xgrid = np.arange(master_dem.nx)
    ygrid = np.arange(master_dem.ny)


    diff_before = master_dem.r-dem2coreg
    if args.plot==True:
      maxval = np.percentile(np.abs(diff_before[np.isfinite(diff_before)]),99)
      pl.imshow(diff_before,vmin=-maxval,vmax=maxval)
      cb=pl.colorbar()
      cb.set_label('Elevation difference (m)')
      pl.show()

    #Create spline functions
    f = RectBivariateSpline(ygrid,xgrid, slave_filled,kx=1,ky=1)
    f2 = RectBivariateSpline(ygrid,xgrid, nanval,kx=1,ky=1)
    xoff, yoff = 0,0 
    
    for i in xrange(args.niter):

        #compute offset
        east, north, c = coregistration(master_dem.r,dem2coreg,args.plot)
        print "Offset in pixels : (%f,%f)" %(east,north)
        xoff+=north
        yoff+=east
    
        #resample slave DEM in the new grid
        znew = f(ygrid-yoff,xgrid-xoff)
        nanval_new = f2(ygrid-yoff,xgrid-xoff)
        
        #remove filled values that have been interpolated
        znew[nanval_new!=0] = np.nan

        #Remove bias
        diff = znew-master_dem.r
    #    bias = np.median(diff[np.isfinite(diff)])
        bias = np.nanmean(diff)
        znew-= bias

        dem2coreg = znew    


    #Display results
    if args.plot==True:
        diff_after = master_dem.r - znew

        pl.figure('before')
        pl.imshow(diff_before,vmin=-maxval,vmax=maxval)
        cb = pl.colorbar()
        cb.set_label('DEM difference (m)')
        pl.figure('after')
        pl.imshow(diff_after,vmin=-maxval,vmax=maxval)
        cb = pl.colorbar()
        cb.set_label('DEM difference (m)')
        #pl.show()

        pl.figure()
        pl.hist(diff_after[np.isfinite(diff_after)],bins=np.linspace(-maxval,maxval,50))
        pl.xlabel('DEM difference (m)')
        pl.show()

    #Save to output file
    dtype = master_dem.ds.GetRasterBand(1).DataType
    raster.simple_write_geotiff(args.outfile, znew, master_dem.ds.GetGeoTransform(), wkt=master_dem.srs.ExportToWkt(),dtype=dtype)
