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
import os, sys
import gdal
from glob import glob

#Personal libraries
import georaster as raster
from geoutils.demraster import DEMRaster
from geoutils import geometry as geo
from geoutils import geovector as vect

#Disable warnings
import warnings
warnings.filterwarnings("ignore")


def grad2d(dem):
  '''
  Calculate the slope and gradient of a DEM
  '''

#  I = ndimage.gaussian_filter(dem,0.333)
  g2, g1 = np.gradient(dem) # in Python, x and y axis reversed

  slope_pix = np.sqrt(g1**2 + g2**2)
  aspect = np.arctan2(-g1,g2)    #aspect=0 when slope facing north
  aspect = aspect+np.pi

  return slope_pix,aspect


def horizontal_shift(dh,slope,aspect,plot=False, min_count=30):
    """
    Compute the horizontal shift between 2 DEMs using the method presented in Nuth & Kaab 2011
    Inputs :
    - dh : array, elevation difference master_dem - slave_dem
    - slope/aspect : array, slope and aspect for the same locations as the dh
    - plot : bool, set to True to display the slope/aspect graph
    - min_count : int, minimum number of points in aspect bins for it to be considered valid (Default is 30)
    Returns :
    - east, north, c : f, estimated easting and northing of the shift, c is not used here but is related to the vertical shift
    """

    # function to be correlated according to Nuth & Kaab
    x = aspect
    y = dh/slope

    #remove non-finite values
    xf = x[(np.isfinite(x)) & (np.isfinite(y))]
    yf = y[(np.isfinite(x)) & (np.isfinite(y))]

    #remove outliers
    p1 = np.percentile(yf,1)
    p99 = np.percentile(yf,99)
    valids = np.where((p1<yf) & (yf<p99) & (np.abs(yf)<200))
    xf = xf[valids]
    yf = yf[valids]

    # compute median value for different aspect slices
    slice_bounds = np.arange(0,2*np.pi,np.pi/36)
    y_med=np.zeros([len(slice_bounds)])
    count=np.zeros([len(slice_bounds)])
    j=0
    for i in slice_bounds:
        yf_slice = yf[(i<xf) & (xf<i+np.pi/36)] # extract y values in bin
        y_med[j] = np.median(yf_slice) # calculate median
        count[j] = len(yf_slice)  # count number of elements
        j=j+1

    # Filter out bins with counts below threshold
    yobs = y_med[count>min_count]
    xs = slice_bounds[count>min_count]

    # Check that enough observations exist
    # Should do something better, like calculating distribution of aspects
    if len(xs)<10:
      print("ERROR: Less than 10 different aspect exist (out of 72), increase min_count or extend the DEMs")
      if plot==True:
        ax1 = pl.subplot(111)
        p1 = ax1.plot(xs,yobs,'b.',label='Obs.',ms=12)
        pl.xlabel('Terrain aspect (rad)')
        pl.ylabel(r'dh/tan($\alpha$)')
        ax2 = ax1.twinx()
        p3 = ax2.bar(slice_bounds,count,width=np.pi/36,facecolor='gray',alpha=0.2, label='count')
        ax2.plot(xs,(min_count,)*len(xs),'k-',alpha=0.2)
        ax2.set_ylabel('count')
        ax2.set_xlim([np.min(slice_bounds),np.max(slice_bounds)])
        ax1.legend(loc='best',numpoints=1)
        pl.show()
        
      sys.exit(1)
      
    #First guess
    p0 = (3*np.std(yobs)/(2**0.5),0,np.mean(yobs))

    #Least square fit   
    def peval(x,p):
        return p[0]*np.cos(p[1]-x) + p[2]

    def residuals(p,y,x):
        err = peval(x,p)-y
        return err

    plsq = leastsq(residuals, p0, args = (yobs,xs),full_output = 1)
    # plsq = leastsq(residuals, p0, args = (yf,xf),full_output = 1)
    yfit = peval(slice_bounds,plsq[0])

    #plotting results
    if plot==True:
        ax1 = pl.subplot(111)
        p1 = ax1.plot(xs,yobs,'b.',label='Obs.',ms=12)
        p2 = ax1.plot(slice_bounds,yfit,'k-',label='Fit',lw=2)
        pl.xlabel('Terrain aspect (rad)')
        pl.ylabel(r'dh/tan($\alpha$)')
        ax2 = ax1.twinx()
        p3 = ax2.bar(slice_bounds,count,width=np.pi/36,facecolor='gray',alpha=0.2, label='count')
        ax2.plot(slice_bounds,(min_count,)*len(slice_bounds),'k-',alpha=0.2)
        ax2.set_ylabel('count')
        ax2.set_xlim([np.min(slice_bounds),np.max(slice_bounds)])
        ax1.legend(loc='best',numpoints=1)
        pl.show()

    a,b,c = plsq[0]
    east = a*np.sin(b)     #with b=0 when north (origin=y-axis)
    north = a*np.cos(b)

    return east, north, c

def deramping(diff,X,Y,d=1,plot=False):
  """
  Estimate a ramp (tilt) in elevation difference between two DEMs.
  Inputs :
  - diff : array, elevation difference between the two DEMs
  - X, Y : arrays, X, Y position of the elevation difference in any system
  - d : degree of the polynomial to remove (deg 1: a0 + a1*X + a2*Y; deg 2: a0 +a1*X +a2*Y + a3*X**2 +a4*X*Y + a5*Y**2)
  - plot : i f set to True, plots are displayed
  Returns :
  - ramp : a function that defines the estimated ramp, if two arguments X,Y are passed, return the value of the ramp at each location
  """

  #filter outliers
  # med = np.median(diff[np.isfinite(diff)])
  # mad=1.4826*np.median(np.abs(diff[np.isfinite(diff)]-med))
  # diff[np.abs(diff)>3*mad] = np.nan

  #Least square fit
  def poly2D(X,Y,p,d):
    """
    Calculate the elements of a 2D polynomial of degree d, with coefficients p. The index of p increases with monomial degree, then with power of X increasing and power of Y decreasing. E.g. for degree 2:
    p[0] + p[1] * X + p[2] * Y + p[3] * X**2 +p[4] * X*Y +p[5] * Y**2
    p must be of size (d+1)*(d+2)/2
    """
    psize = (d+1)*(d+2)/2
    if len(p)!=psize:
      print("ERROR: p must be of size (d+1)*(d+2)/2 = %i" %psize)
      return None
    else:
      # k is degree considered, j is power of Y varying from 0 to k, k-j is power of X varying from k to 0, k*(k+1)/2 + j is index of the coefficient
      return np.sum([p[k*(k+1)//2+j]*X**(k-j)*Y**j for k in range(d+1) for j in range(k+1)],axis=0)
    

  def residuals(p,z,X,Y,d):
    err = poly2D(X,Y,p,d)-z
    err = err[np.isfinite(err)]
    return err

  z = diff[np.isfinite(diff)]
  x = X[np.isfinite(diff)]
  y = Y[np.isfinite(diff)]

  # reduce number of elements for speed
  if len(x)>5e5:
    inds = np.random.randint(0,len(x)-1,int(5e5))
    x = x[inds]
    y = y[inds]
    z = z[inds]

  p0 = np.zeros((d+1)*(d+2)//2)   # inital guess for the polynomial coeff
  plsq = leastsq(residuals, p0, args = (z,x,y,d),full_output = 1)
  zfit = poly2D(X,Y,plsq[0],d)
  
  if plot==True:
    vmax1 = np.nanmax(np.abs(diff))
    vmax2 = np.nanmax(np.abs(diff-zfit))
    pl.figure('before')
    pl.imshow(diff,vmin=-vmax1,vmax=vmax1,cmap='RdYlBu')
    pl.colorbar()
    pl.figure('after')
    pl.imshow(diff-zfit,vmin=-vmax2,vmax=vmax2,cmap='RdYlBu')
    pl.colorbar()
    pl.figure('ramp')
    pl.imshow(zfit,cmap='RdYlBu',vmin=-vmax1,vmax=vmax1)
    pl.colorbar()
    pl.show()

  def ramp(X,Y):
    return poly2D(X,Y,plsq[0],d)

  return ramp


def coreg_with_master_dem(args):
    """
    Coregistration with the use of a master DEM
    """

    ## Read DEMs ##
    # master
    master_dem = DEMRaster(args.master_dem)
    master_dem.r = np.float32(master_dem.r)
    if args.nodata1!='none':
        nodata1 = float(args.nodata1)
    else:
        band=master_dem.ds.GetRasterBand(1)
        nodata1 = band.GetNoDataValue()

    # slave
    slave_dem = raster.SingleBandRaster(args.slave_dem)
    slave_dem.r = np.float32(slave_dem.r)
    if args.nodata2!='none':
      nodata2 = float(args.nodata2)
    else:
      band=slave_dem.ds.GetRasterBand(1)
      nodata2 = band.GetNoDataValue()

    # master, read intersection only
    # extent = slave_dem.intersection(args.master_dem)
    # master_dem = DEMRaster(args.master_dem,load_data=extent, latlon=False)
    # master_dem.r = np.float32(master_dem.r)
    # if args.nodata1!='none':
    #     master_dem.r[master_dem.r==float(args.nodata1)] = np.nan
    # else:
    #     band=master_dem.ds.GetRasterBand(1)
    #     nodata1 = band.GetNoDataValue()
    #     master_dem.r[master_dem.r==nodata1] = np.nan

    ## reproject slave DEM into the master DEM spatial reference system ##
    if master_dem.r.shape!=slave_dem.r.shape:
      if args.grid=='master':
        
        print("Reproject slave DEM")
        dem2coreg = slave_dem.reproject(master_dem.srs, master_dem.nx, master_dem.ny, master_dem.extent[0], master_dem.extent[3], master_dem.xres, master_dem.yres, dtype=6, nodata=nodata2, interp_type=gdal.GRA_Bilinear,progress=True).r
        gt = (master_dem.extent[0], master_dem.xres, 0.0, master_dem.extent[3], 0.0, master_dem.yres)  # Save GeoTransform for saving
        
      elif args.grid=='slave':
        
        print("Reproject master DEM")
        master_dem = master_dem.reproject(slave_dem.srs, slave_dem.nx, slave_dem.ny, slave_dem.extent[0], slave_dem.extent[3], slave_dem.xres, slave_dem.yres, dtype=6, nodata=nodata1, interp_type=gdal.GRA_Bilinear,progress=True)
        #master_dem = DEMRaster(master_dem.ds)  # convert it to a DEMRaster object for later use of specific functions
        dem2coreg=slave_dem.r
        gt = (slave_dem.extent[0], slave_dem.xres, 0.0, slave_dem.extent[3], 0.0, slave_dem.yres)  # Save GeoTransform for saving
        
      else:
        sys.exit("Error : grid must be 'master' or 'slave'")
        
    else:
      dem2coreg = slave_dem.r
      gt = (master_dem.extent[0], master_dem.xres, 0.0, master_dem.extent[3], 0.0, master_dem.yres)  # Save GeoTransform for saving

    if nodata1 is not None: master_dem.r[master_dem.r==nodata1] = np.nan
    if nodata2 is not None: dem2coreg[dem2coreg==nodata2] = np.nan

    ## mask points ##
    if args.maskfile!='none':
        mask = raster.SingleBandRaster(args.maskfile)
        if master_dem.r.shape!=mask.r.shape:
          print("Reproject mask")
          mask = mask.reproject(master_dem.srs, master_dem.nx, master_dem.ny, master_dem.extent[0], master_dem.extent[3], master_dem.xres, master_dem.yres, dtype=6, interp_type=gdal.GRA_NearestNeighbour,progress=True)  # nearest neighbor interpolation

        master_dem.r[mask.r>0] = np.nan

    if args.shp!='none':
      outlines = vect.SingleLayerVector(args.shp)
      mask = outlines.create_mask(master_dem)
      master_dem.r[mask>0] = np.nan

    ## filter outliers ##
    if args.resmax!='none':
      master_dem.r[np.abs(master_dem.r-dem2coreg)>float(args.resmax)] = np.nan

    ## Set master DEM grid for later resampling ##
    xgrid = np.arange(master_dem.nx)
    ygrid = np.arange(master_dem.ny)
    X, Y = master_dem.coordinates()


    diff_before = dem2coreg-master_dem.r


    ## Print out some statistics
    median = np.median(diff_before[np.isfinite(diff_before)])
    NMAD_old = 1.4826*np.median(np.abs(diff_before[np.isfinite(diff_before)]-median))
    print("Statistics on initial dh")
    print("Median : %f, NMAD : %f" %(median,NMAD_old))

    ## Display
    if args.plot==True:
      maxval = 3*NMAD_old #np.percentile(np.abs(diff_before[np.isfinite(diff_before)]),90)
      pl.imshow(diff_before,vmin=-maxval,vmax=maxval,cmap='RdYlBu')
      cb=pl.colorbar()
      cb.set_label('Elevation difference (m)')
      pl.show()

    
    ## fill NaN values for interpolation ##
    nanval = np.isnan(dem2coreg)
    slave_filled = np.where(np.isnan(dem2coreg),-9999,dem2coreg)
    
    ## Create spline function ##
    f = RectBivariateSpline(ygrid,xgrid, slave_filled,kx=1,ky=1)
    f2 = RectBivariateSpline(ygrid,xgrid, nanval,kx=1,ky=1)
    xoff, yoff = 0,0 
    
    ## compute terrain aspect/slope ##
    slope, aspect = grad2d(master_dem.r)


    ## Iterations to estimate DEMs shift
    print("Iteratively estimate DEMs shift")

    for i in range(args.niter):

	# remove bias
        dem2coreg-=median

        #Elevation difference
        dh = master_dem.r - dem2coreg

        #compute offset
        east, north, c = horizontal_shift(dh,slope,aspect,args.plot,args.min_count)
        print("#%i - Offset in pixels : (%f,%f)" %(i+1,east,north))
        xoff+=east
        yoff+=north
    
        #resample slave DEM in the new grid
        znew = f(ygrid-yoff,xgrid+xoff)    #postive y shift moves south
        nanval_new = f2(ygrid-yoff,xgrid+xoff)
        
        #remove filled values that have been interpolated
        znew[nanval_new!=0] = np.nan

	# update DEM
        dem2coreg = znew    

	# print some statistics
        diff = dem2coreg-master_dem.r
        diff = diff[np.isfinite(diff)]
        NMAD_new = 1.4826*np.median(np.abs(diff-np.median(diff)))
        median = np.median(diff)

        print("Median : %.2f, NMAD = %.2f, Gain : %.2f%%" %(median,NMAD_new,(NMAD_new-NMAD_old)/NMAD_old*100))
        NMAD_old = NMAD_new

    print("Final Offset in pixels (east, north) : (%f,%f)" %(xoff,yoff))
    
    if args.save==True:
      fname, ext = os.path.splitext(args.outfile)
      fname+='_shift.txt'
      f = open(fname,'w')
      f.write("Final Offset in pixels (east, north) : (%f,%f)" %(xoff,yoff))
      f.write("Final NMAD : %f" %NMAD_new)
      f.close()
      print("Offset saved in %s" %fname)
      
    ### Deramping ###
    print("Deramping")
    diff = dem2coreg-master_dem.r
    
    # remove points above altitude threshold (snow covered areas) 
    if args.zmax!='none':
      diff[master_dem.r>int(args.zmax)] = np.nan

    # remove points below altitude threshold (e.g sea ice)
    if args.zmin!='none':
      diff[master_dem.r<int(args.zmin)] = np.nan

    # remove points with slope higher than 40° that are more error-prone
    # removed until issue with Nan is fixed, see previous commit
    #slope, aspect = master_dem.compute_slope()
    #diff[slope>=40*np.pi/180] = np.nan
    #diff[np.isnan(slope)] = np.nan

    # remove outliers
    med = np.median(diff[np.isfinite(diff)])
    mad=1.4826*np.median(np.abs(diff[np.isfinite(diff)]-med))
    diff[np.abs(diff-med)>3*mad] = np.nan

    # estimate a ramp and remove it
    ramp = deramping(diff,X,Y,d=args.degree,plot=args.plot)
    dem2coreg-=ramp(X,Y)

    # save to output file
    if args.save==True:
      fname, ext = os.path.splitext(args.outfile)
      fname+='_ramp.TIF'
      #fname = WD+'/ramp.out'
      raster.simple_write_geotiff(fname, ramp(X,Y), gt, wkt=master_dem.srs.ExportToWkt(),dtype=gdal.GDT_Float32)
      #ramp(X,Y).tofile(fname)
      print("Ramp saved in %s" %fname)
      
    # print some statistics
    diff = dem2coreg-master_dem.r
    diff = diff[np.isfinite(diff)]
    median = np.median(diff)
    NMAD = 1.4826*np.median(np.abs(diff-median))
    print("Final DEM")
    print("Median : %.2f, NMAD = %.2f" %(median,NMAD))


    #Display results
    if args.plot==True:
        diff_after = dem2coreg - master_dem.r

        pl.figure('before')
        pl.imshow(diff_before,vmin=-maxval,vmax=maxval,cmap='RdYlBu')
        cb = pl.colorbar()
        cb.set_label('DEM difference (m)')
        pl.figure('after')
        pl.imshow(diff_after,vmin=-maxval,vmax=maxval,cmap='RdYlBu')
        cb = pl.colorbar()
        cb.set_label('DEM difference (m)')
        #pl.show()

        pl.figure()
        pl.hist(diff_after[np.isfinite(diff_after)],bins=np.linspace(-maxval,maxval,50))
        pl.xlabel('DEM difference (m)')
        pl.show()

    #Save to output file
    #dtype = master_dem.ds.GetRasterBand(1).DataType
    raster.simple_write_geotiff(args.outfile, dem2coreg, gt, wkt=master_dem.srs.ExportToWkt(),dtype=gdal.GDT_Float32)


def read_icesat_elev(is_files,RoI):
    """
    Read IceSAT elevation in a series of files, within a region of interest RoI
    Inputs :
    - is_files : str, IceSAT file name or regular expression to different files
    - RoI : str, region of interest, defined as a polygon ((x1,y1),(x2,y2),...)
    """

    import h5py

    ## Read Icesat data
    print("read Icesat data")
    filelist = glob(is_files)
    filelist.sort()

    all_lats = []
    all_lons = []
    all_elev = []

    for f in filelist:
        #print(f)
        ds=h5py.File(f,'r')
        lons = np.array(ds['Data_40HZ/Geolocation/d_lon'])
        lats = np.array(ds['Data_40HZ/Geolocation/d_lat'])

        inside=geo.points_inside_polygon(lons,lats,RoI)
        if len(inside[inside==True]):
            lons = lons[inside==True]
            lats = lats[inside==True]
            elev = np.array(ds['Data_40HZ/Elevation_Surfaces/d_elev'][inside==True])
            all_lons.extend(lons)
            all_lats.extend(lats)
            all_elev.extend(elev)

    all_lons = np.array(all_lons)
    all_lats = np.array(all_lats)
    all_elev = np.array(all_elev)
    
    return all_lons, all_lats, all_elev


def coreg_with_IceSAT(args):
    """
    Coregistration with the use of IceSAT data
    """
    is_files = args.master_dem

    ## Read slave DEM ##
    slave_dem = DEMRaster(args.slave_dem)
    slave_dem.r = np.float32(slave_dem.r)
    if args.nodata2!='none':
      nodata = float(args.nodata2)
    else:
      band=slave_dem.ds.GetRasterBand(1)
      nodata = band.GetNoDataValue()

    # save original data for later resampling
    dem2coreg = slave_dem
    dem2coreg_save = np.copy(dem2coreg.r)

    # compute DEM extent
    lonmin,lonmax,latmin,latmax=dem2coreg.get_extent_latlon()
    lonmin+=360
    lonmax+=360
    RoI = ((lonmin,latmin),(lonmin,latmax),(lonmax,latmax),(lonmax,latmin),(lonmin,latmin))

    ## mask points ##
    mask = raster.SingleBandRaster(args.maskfile)
    if dem2coreg.r.shape!=mask.r.shape:
          print("Reproject mask")
          mask = mask.reproject(dem2coreg.srs, dem2coreg.nx, dem2coreg.ny, dem2coreg.extent[0], dem2coreg.extent[3], dem2coreg.xres, dem2coreg.yres, dtype=6, nodata=nodata, interp_type=1,progress=True)

    dem2coreg.r[mask.r>0] = np.nan

    ## read Icesat data with DEM extent ##
    all_lons, all_lats, all_elev = read_icesat_elev(is_files,RoI)

    ## compare slave DEM and Icesat ##
    slave_elev=dem2coreg.interp(all_lons,all_lats,latlon=True)
    dh = all_elev - slave_elev
    dh[slave_elev==0] = np.nan
    xx, yy = dem2coreg.proj(all_lons,all_lats)

    if args.plot==True:
        pl.title('Icesat - slave DEM elev')
        rgb=dem2coreg.shaded_relief(downsampl=5)
        pl.scatter(xx,yy,c=dh,edgecolor='none',vmin=-20,vmax=20)
        cb=pl.colorbar()
        cb.set_label('Elevation difference (m)')
        pl.show()

    ## compute slave DEM slope at Icesat points ##
    print("Compute slope and aspect")
    g2, g1 = np.gradient(dem2coreg.r)
    distx = np.abs(dem2coreg.xres)
    disty = np.abs(dem2coreg.yres)
    slope_pix = np.sqrt((g1/distx)**2 + (g2/disty)**2)
    aspect = np.arctan2(-g1,g2)
    aspect = aspect+np.pi

    slope_ds = raster.simple_write_geotiff('none', slope_pix, dem2coreg.ds.GetGeoTransform(), wkt=dem2coreg.srs.ExportToWkt(), dtype=6)
    slope_raster = raster.SingleBandRaster(slope_ds)
    slope_at_IS = slope_raster.interp(all_lons,all_lats,latlon=True)

    aspect_ds = raster.simple_write_geotiff('none', aspect, dem2coreg.ds.GetGeoTransform(), wkt=dem2coreg.srs.ExportToWkt(), dtype=6)
    aspect_raster = raster.SingleBandRaster(aspect_ds)
    aspect_at_IS = aspect_raster.interp(all_lons,all_lats,latlon=True)

    # slave DEM grid
    xgrid = np.arange(dem2coreg.nx)
    ygrid = np.arange(dem2coreg.ny)
    X, Y = dem2coreg.coordinates()

    ## Print out some statistics
    median = np.median(dh[np.isfinite(dh)])
    NMAD_old = 1.4826*np.median(np.abs(dh[np.isfinite(dh)]-median))
    print("Statistics on initial dh")
    print("Median : %f, NMAD : %f" %(median,NMAD_old))



    ## Iterations to estimate DEMs shift
    print("Iteratively estimate DEMs shift")

    slave_elev=dem2coreg.interp(all_lons,all_lats,latlon=True)
    dh = all_elev - slave_elev
    dh[slave_elev==0] = np.nan
    xoff, yoff = 0,0 

    for i in range(args.niter):
	
        # compute aspect/dh relationship
        east, north, c = horizontal_shift(dh,slope_at_IS,aspect_at_IS,plot=args.plot, min_count=args.min_count)
        print("#%i - Offset in pixels : (%f,%f)" %(i+1,east,north))
        xoff+=east
        yoff+=north

        #Update elevation difference
        slave_elev=dem2coreg.interp(xx+xoff,yy+yoff)
        dh = all_elev-slave_elev
        dh[slave_elev==0] = np.nan

	# print some statistics
        median = np.median(dh[np.isfinite(dh)])
        NMAD_new = 1.4826*np.median(np.abs(dh[np.isfinite(dh)]-median))

        print("Median : %.2f, NMAD = %.2f, Gain : %.2f%%" %(median,NMAD_new,(NMAD_new-NMAD_old)/NMAD_old*100))
        NMAD_old = NMAD_new

    print("Final Offset in pixels (east, north) : (%f,%f)" %(xoff,yoff))

    if args.save==True:
      fname, ext = os.path.splitext(args.outfile)
      fname+='_shift.txt'
      f = open(fname,'w')
      f.write("Final Offset in pixels (east, north) : (%f,%f)" %(xoff,yoff))
      f.write("Final NMAD : %f" %NMAD_new)
      f.close()
      print("Offset saved in %s" %fname)
      

    ### Deramping ###
    print("Deramping")
    
    # remove points above altitude threshold (snow covered areas)
    #if args.zmax!='none':
    #  dh[master_dem.r>int(args.zmax)] = np.nan

    # remove points below altitude threshold (e.g sea ice)
    zmin=40
    if zmin!='none':
      dh[slave_elev<int(zmin)] = np.nan

    # remove points with slope higher than 20° that are more error-prone
#    slope, aspect = dem2coreg.compute_slope()
#    dh[slope>=20*np.pi/180] = np.nan
#    dh[np.isnan(slope)] = np.nan

    # remove outliers
    med = np.median(dh[np.isfinite(dh)])
    mad=1.4826*np.median(np.abs(dh[np.isfinite(dh)]-med))
    dh[np.abs(dh-med)>3*mad] = np.nan

    # estimate a ramp and remove it
    ramp = deramping(dh,xx,yy,d=args.degree,plot=False)

    # compute stats of deramped dh
    tmp = dh-ramp(X,Y) ;
    median = np.median(tmp)
    NMAD_new = 1.4826*np.median(np.abs(tmp[np.isfinite(tmp)]-median))
    
    if args.save==True:
          fname, ext = os.path.splitext(args.outfile)
          fname+='_shift.txt'
          f = open(fname,'a')
          f.write("Median after deramping : (%f)" %(median))
          f.write("NMAD after deramping : %f" %NMAD_new)
          f.close()
          print("Post-deramping stats saved in %s" %fname)
        
    # save to output file
    if args.save==True:
      fname, ext = os.path.splitext(args.outfile)
      fname+='_ramp.TIF'
      #fname = WD+'/ramp.out'
      raster.simple_write_geotiff(fname, ramp(X,Y), dem2coreg.ds.GetGeoTransform(), wkt=dem2coreg.srs.ExportToWkt(),dtype=gdal.GDT_Float32)
      #ramp(X,Y).tofile(fname)
      print("Ramp saved in %s" %fname)
      
    if args.plot==True:
        pl.figure('ramp')
        pl.scatter(xx,yy,c=ramp(xx,yy),edgecolor='none')
        pl.colorbar()
        pl.figure('before')
        pl.imshow(rgb,extent=dem2coreg.extent,interpolation='bilinear')
        pl.scatter(xx,yy,c=dh,edgecolor='none',vmin=-10,vmax=10)
        pl.colorbar()
        pl.figure('after')
        pl.imshow(rgb,extent=dem2coreg.extent,interpolation='bilinear')
        pl.scatter(xx,yy,c=dh-ramp(xx,yy),edgecolor='none',vmin=-10,vmax=10)
        pl.colorbar()
        pl.show()

    
    ### Interpolate the slave DEM to the new grid ###

    print("Interpolate DEM to new grid")

    # fill NaN values for interpolation
    nanval = np.isnan(dem2coreg_save)
    slave_filled = np.where(np.isnan(dem2coreg_save),-9999,dem2coreg_save)

    # Create spline function
    f = RectBivariateSpline(ygrid,xgrid, slave_filled,kx=1,ky=1)
    f2 = RectBivariateSpline(ygrid,xgrid, nanval,kx=1,ky=1)

    # resample slave DEM in the new grid
    znew = f(ygrid-yoff,xgrid+xoff)    #postive y shift moves south
    nanval_new = f2(ygrid-yoff,xgrid+xoff)
        
    #remove filled values that have been interpolated
    znew[nanval_new!=0] = np.nan

    # update DEM
    dem2coreg_save = znew    


    ### Remove ramp ###
    dem2coreg_save-=ramp(X,Y)


    ### Save to output file ###
    raster.simple_write_geotiff(args.outfile,dem2coreg_save,dem2coreg.ds.GetGeoTransform(),wkt=dem2coreg.srs.ExportToWkt(),dtype=gdal.GDT_Float32)


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
    parser.add_argument('-m', dest='maskfile', type=str, default='none', help='str, path to a mask of same size as the master DEM, to filter out non stable areas such as glaciers. Points with mask>0 are masked.  (default is none)')
    parser.add_argument('-shp', dest='shp', type=str, default='none', help='str, path to a shapefile containing outlines of objects to mask, such as RGI shapefiles (default is none)')
    parser.add_argument('-n1', dest='nodata1', type=str, default='none', help='int, no data value for master DEM if not specified in the raster file (default read in the raster file)')
    parser.add_argument('-n2', dest='nodata2', type=str, default='none', help='int, no data value for slave DEM if not specified in the raster file (default read in the raster file)')
    parser.add_argument('-min_count', dest='min_count', type=int, default=30, help='int, minimum number of points in each aspect bin to be considered valid (default is 30)')
    parser.add_argument('-zmax', dest='zmax', type=str, default='none', help='float, points with altitude above zmax are masked during the vertical alignment, e.g snow covered areas (default none)')
    parser.add_argument('-zmin', dest='zmin', type=str, default='none', help='float, points with altitude below zmin are masked during the vertical alignment, e.g points on sea (default none)')
    parser.add_argument('-resmax', dest='resmax', type=str, default='none', help='float, maximum value of the residuals, points where |dh|>resmax are considered as outliers and removed (default none)')
    parser.add_argument('-deg', dest='degree', type=int, default=1, help='int, degree of the polynomial to be fit to residuals and removed (default is 1=tilt)')
    parser.add_argument('-grid', dest='grid', type=str, default='master', help="'master' or 'slave' : common grid to use for the DEMs (default is master DEM grid)")
    parser.add_argument('-save', dest='save', help='Save horizontal offset as a text file and ramp as a GTiff file',action='store_true')
    parser.add_argument('-IS', dest='IS', help='Master DEM are IceSAT data instead of a raster DEM. master_dem must then be a string to the file name or regular expression to several files (use quotes)',action='store_true')


    args = parser.parse_args()

    if args.IS==False:
      coreg_with_master_dem(args)
    else:
      coreg_with_IceSAT(args)
    
