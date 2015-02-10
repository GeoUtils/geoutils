#!/usr/bin/env python
#coding=utf-8

from osgeo import ogr, osr, gdal
import matplotlib.pyplot as pl
import georaster as raster
import numpy as np
import numbers
from copy import deepcopy

"""
geovector.py

Classes to simplify reading geographic vector data and associated 
georeferencing information.

There are two classes available:

    SingleLayerVector    --  For loading datasets with only a single layer of 
                            data
    MultiLayerVector     --  For loading datasets that contain multiple layers 
                            of data (TO IMPLEMENT)

Each class works as a wrapper to the GDAL API. A vector dataset can be loaded
into a class without needing to load the actual data as well, which is useful
for querying geo-referencing information without memory overheads.

Both classes provide comprehensive georeferencing information access via
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



Examples
--------

...


Created on Mon Feb 09 2015

@author: Amaury Dehecq (amaury.dehecq@univ-savoie.fr)
"""


# By default, GDAL/OGR does not raise exceptions - enable them
# See http://trac.osgeo.org/gdal/wiki/PythonGotchas
gdal.UseExceptions()
ogr.UseExceptions()


class Object(object):
    pass


class SingleLayerVector:

    def __init__(self,ds_filename,load_data=True,latlon=True,band=1):
        """ Construct object from vector file with a single layer. 
        
        Parameters:
            ds_filename : filename of the dataset to import
            
        """
        self._load_ds(ds_filename)     



    def __del__(self):
        """ Close the gdal link to the dataset on object destruction """
        self.ds = None



    def _load_ds(self,ds_filename):
        """ Load link to data file and set up georeferencing """
        self.ds_file = ds_filename
        self.ds = ogr.Open(ds_filename)

        # Check to see if shapefile is found.
        if self.ds is None:
            print 'Could not open %s' % (ds_filename)

        self.lyr = self.ds.GetLayer()
#        self.ds.lyr = []
#        for i in xrange(self.ds.GetLayerCount()):
#            self.ds.lyr.append(self.ds.GetLayerByIndex)

        self.srs = self.lyr.GetSpatialRef()
        self.extent = self.lyr.GetExtent()

        layerDefinition = self.lyr.GetLayerDefn()
        
        fields = []
        dtypes = []
        for i in range(layerDefinition.GetFieldCount()):
            fields.append(layerDefinition.GetFieldDefn(i).GetName())
            fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
            dtypes.append(layerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode))
        """
        feat = self.lyr.GetFeature(0)
        fields = feat.keys()
        
        dtypes = []
        for i in xrange(len(fields)):
            dt = feat.GetFieldType(fields[i])
            dtypes.append(gdal.GetDataTypeName(dt))
        """
        self.fields = fields
        self.dtypes = dtypes
        
    def read(self):
        """
        Read features defined in the vector file.
        Features filtered before are not read.
        """

        features = []
        for i in xrange(self.lyr.GetFeatureCount()):
            features.append(self.lyr.GetNextFeature())
        self.features = np.array(features)


    def plot(self,indices='all'):
        """
        Plot the geometries defined in the vector file
        Inputs :
        - indices : indices of the features to plot (Default is 'all')
        """

        if not hasattr(self,'features'):
            print "You must run self.read() first"
            return 0

        if indices=='all':
            indices = range(len(self.features))
        elif isinstance(indices,numbers.Number):
            indices = [indices,] #create list if only one value

        for feat in self.features[indices]:
            sh = Shape(feat)
            sh.plot()


    def fill(self,indices='all'):
        """
        Display the filled geometries defined in the vector file
        Inputs :
        - indices : indices of the features to plot (Default is 'all')
        """

        if not hasattr(self,'features'):
            print "You must run self.read() first"
            return 0

        if indices=='all':
            indices = range(len(self.features))
        elif isinstance(indices,numbers.Number):
            indices = [indices,] #create list if only one value

        for feat in self.features[indices]:
            sh = Shape(feat)
            sh.fill()


    def crop(self,xmin,xmax,ymin,ymax):
        """
        Filter all features that are outside the rectangle defined by the corners xmin, xmax, ymin, ymax.
        Coordinates are in the same Reference System as the vector file.
        """
        # Create a ring linking the 4 corners
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(xmin, ymin)
        ring.AddPoint(xmin, ymax)
        ring.AddPoint(xmax, ymax)
        ring.AddPoint(xmax, ymin)
        ring.AddPoint(xmin, ymin)

        # Create polygon
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        #other method
#        wkt = "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))" %(xmin,ymin,xmin,ymax,xmax,ymax,xmax,ymin,xmin,ymin)
#        poly=ogr.CreateGeometryFromWkt(wkt)
        
        #Filter features outside the polygon
        self.lyr.SetSpatialFilter(poly)

    
    def crop2raster(self,rasterfile):
        """
        Filter all features that are outside the raster extent.
        """

        # Read raster extent
        img = raster.SingleBandRaster(rasterfile,load_data=False)
        xmin, xmax, ymin, ymax = img.extent

        # Reproject extent to vector geometry
        sourceSR = img.srs
        targetSR = self.srs
        transform = osr.CoordinateTransformation(sourceSR,targetSR)
        left, bottom = transform.TransformPoint(xmin,ymin)[:2]
        right, up = transform.TransformPoint(xmax,ymax)[:2]
        
        # Crop
        self.crop(left,right,bottom,up)

    def zonal_statistics(self,rasterfile,indices='all',nodata=None):
        """
        Compute statistics of the data in rasterfile for each feature in self.
        Inputs :
        - indices : indices of the features to compute (Default is 'all')
        """

        #Dimension vector file to raster size
#        self.crop2raster(rasterfile)
#        self.read()

        # Read raster coordinates
        img = raster.SingleBandRaster(rasterfile)
        X, Y = img.coordinates()
#        left, right, bottom, up = img.extent
#        xsize, ysize = img.get_pixel_size()
#        X, Y = np.meshgrid(np.arange(left,right,abs(xsize)),np.arange(bottom,up,abs(ysize)))

        # Reproject vector geometry to same projection as raster
        sourceSR = self.srs
        targetSR = img.srs
        coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)

        if indices=='all':
            indices = range(len(self.features))
        elif isinstance(indices,numbers.Number):
            indices = [indices,] #create list if only one value

        median = []
        std = []
        count = []
        for feat in self.features[indices]:
            sh = Shape(feat)
            geom = deepcopy(sh.geom)
            geom.Transform(coordTrans)
            xmin, xmax, ymin, ymax = geom.GetEnvelope()

            inds = np.where((X>=xmin) & (X<=xmax) & (Y>=ymin) & (Y<=ymax))

            inside_i, inside_j = [], []
            for i,j in np.transpose(inds):
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(X[i,j], Y[i,j])
                if point.Within(geom):
                    inside_i.append(i)
                    inside_j.append(j)
            
            data = img.r[inside_i,inside_j]
            data = data[~np.isnan(data)]
            if nodata!=None:
                data = data[data!=nodata]

            if len(data)>0:
                median.append(np.median(data))
                std.append(np.std(data))
                count.append(len(data))
            else:
                median.append(np.nan)
                std.append(np.nan)
                count.append(0)

        return np.array(median), np.array(std), np.array(count)


class Shape():

        def __init__(self,feature):
            
            self.feature=feature
            self.geom = feature.GetGeometryRef()
            self.extent = self.geom.GetEnvelope()
            
            self.points = Object()
            self.points.x = []
            self.points.y = []
            if self.geom.GetGeometryCount()>0:
                for i in xrange(self.geom.GetGeometryCount()):
                    poly = self.geom.GetGeometryRef(i)
                    x, y = [], []
                    for j in xrange(poly.GetPointCount()):
                        x.append(poly.GetX(j))
                        y.append(poly.GetY(j))
                 
                    self.points.x.append(x)
                    self.points.y.append(y)

        def plot(self):

            for k in xrange(len(self.points.x)):
                pl.plot(self.points.x[k],self.points.y[k],'b')


        def fill(self):

            if len(self.points.x)>1:
                pl.fill(self.points.x[0],self.points.y[0],'b')

                for k in xrange(1,len(self.points.x)):
                    pl.fill(self.points.x[k],self.points.y[k],'w')
            else:
                pl.fill(self.points.x,self.points.y,'b')
