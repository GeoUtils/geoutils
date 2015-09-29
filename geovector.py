#!/usr/bin/env python
#coding=utf-8

#Python libraries
from osgeo import ogr, osr, gdal
import matplotlib.pyplot as pl
import numpy as np
import numbers
from copy import deepcopy
import mpl_toolkits.basemap.pyproj as pyproj
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib import cm
from matplotlib.collections import PatchCollection
from matplotlib import colors
from scipy import ndimage
import sys
try:
    from skimage import morphology
except ImportError:
    pass

#Personal libraries
import georaster as raster
import geometry as geo

"""
geovector.py

Classes to simplify reading geographic vector data and associated 
georeferencing information.

There are two classes available:

    SingleLayerVector    --  For loading datasets with only a single layer of 
                            data
    MultiLayerVector     --  For loading datasets that contain multiple layers 
                            of data (TO IMPLEMENT)

       The class works as a wrapper to the OGR API. A vector dataset can be loaded
into a class without needing to load the actual data, which is useful
for querying geo-referencing information without memory overheads.

The following attributes are available :

    class.ds    :   the OGR handle to the dataset, which provides access to 
                    all OGR functions, e.g. GetLayer, GetLayerCount...
                        More information on API:
                        http://gdal.org/python/osgeo.ogr.DataSource-class.html

    class.srs   :   an OGR Spatial Reference object representation of the 
                    dataset.
                        More information on API:
                        http://www.gdal.org/classOGRSpatialReference.html

    class.proj  :   a pyproj coordinate conversion function between the 
                    dataset coordinate system and lat/lon.

    class.extent :  tuple of the corners of the dataset in native coordinate
                    system, as (left,right,bottom,top).

    class.extent :  tuple of the corners of the dataset in native coordinate
                    system, as (left,right,bottom,top).

    class.layer : an OGR Layer object, see http://gdal.org/python/osgeo.ogr.Layer-class.html

    class.fields : an OGR field, i.e attributes to each feature. 
                   - fields.name contains the name of the fields
                   - fields.dtype contains the Numpy.dtype of the fields
                   - fields.value contains the value of the fields in form of a dictionary

Additionally, a number of instances are available in the class. 



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

       The class works as a wrapper to the OGR API. A vector dataset can be loaded
into a class without needing to load the actual data, which is useful
for querying geo-referencing information without memory overheads.

The following attributes are available :

    class.ds    :   the OGR handle to the dataset, which provides access to 
                    all OGR functions, e.g. GetLayer, GetLayerCount...
                        More information on API:
                        http://gdal.org/python/osgeo.ogr.DataSource-class.html

    class.srs   :   an OGR Spatial Reference object representation of the 
                    dataset.
                        More information on API:
                        http://www.gdal.org/classOGRSpatialReference.html

    class.proj  :   a pyproj coordinate conversion function between the 
                    dataset coordinate system and lat/lon.

    class.extent :  tuple of the corners of the dataset in native coordinate
                    system, as (left,right,bottom,top).

    class.extent :  tuple of the corners of the dataset in native coordinate
                    system, as (left,right,bottom,top).

    class.layer : an OGR Layer object, see http://gdal.org/python/osgeo.ogr.Layer-class.html

    class.fields : an OGR field, i.e attributes to each feature. 
                   - fields.name contains the name of the fields
                   - fields.dtype contains the Numpy.dtype of the fields
                   - fields.value contains the value of the fields in form of a dictionary

Additionally, a number of instances are available in the class. 
        """

        self._load_ds(ds_filename)     



    def __del__(self):
        """ Close the gdal link to the dataset on object destruction """
        self.ds = None



    def _load_ds(self,ds_filename):
        """ Load link to data file and set up georeferencing """
        if isinstance(ds_filename,str):
            self.ds_file = ds_filename
            self.ds = ogr.Open(ds_filename,0) #read-only
        elif isinstance(ds_filename,ogr.DataSource):
            self.ds = ds_filename
            self.ds_file = ds_filename.GetName()

        # Check to see if shapefile is found.
        if self.ds is None:
            print 'Could not open %s' % (ds_filename)

        #Get layer
        self.layer = self.ds.GetLayer()
        
        # For multiple layers
        # self.ds.layer = []
        # for i in xrange(self.ds.GetLayerCount()):
        # self.ds.layer.append(self.ds.GetLayerByIndex)

        # Get georeferencing infos
        self.srs = self.layer.GetSpatialRef()
        self.extent = self.layer.GetExtent()

        # Get layer fields name and type
        layerDefinition = self.layer.GetLayerDefn()
        
        fields = Object()
        fields.name = []
        fields.dtype = []
        fields._size = []
        for i in range(layerDefinition.GetFieldCount()):
            fields.name.append(layerDefinition.GetFieldDefn(i).GetName())
            fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
            fields.dtype.append(layerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode))
            fields._size.append(layerDefinition.GetFieldDefn(i).GetWidth())

        # Convert types to Numpy dtypes
        # see http://www.gdal.org/ogr__core_8h.html#a787194bea637faf12d61643124a7c9fc
        # IntegerList, RealList, StringList, Integer64List not implemented yet.
        OGR2np = {'Integer':'i4','Real':'d','String':'S','Binary':'S','Date':'S','Time':'S','DateTime':'S','Integer64':'i8'}

        for k in xrange(len(fields.dtype)):
            dtype = OGR2np[fields.dtype[k]]
            if dtype=='S':
                dtype+=str(fields._size[k])
            fields.dtype[k] = dtype

        self.fields = fields
        
        try:
            self.proj = pyproj.Proj(self.srs.ExportToProj4())
        except AttributeError:  #case srs not defined
            self.proj = None


        
    def read(self,subset='all'):
        """
        Load features and fields values defined in the vector file.
        Features filtered before are not read.
        """

        if subset=='all':
            nFeat = self.layer.GetFeatureCount()
        else:
            nFeat = len(subset)
        
        #Initialize dictionnary containing fields data
        self.fields.values = {}
        for k in xrange(len(self.fields.name)):
            f = self.fields.name[k]
            dtype = self.fields.dtype[k]
            self.fields.values[f] = np.empty(nFeat,dtype=dtype)

        if subset!='all':
            if isinstance(subset,numbers.Number):
                subset = [subset,] #create list if only one value

        features = []
        if subset=='all':
            for i in xrange(nFeat):
                #read each feature
                feat = self.layer.GetNextFeature()
                features.append(feat)

                #read each field associated to the feature
                for f in self.fields.name:
                    self.fields.values[f][i] = feat.GetField(f)
        else:
            k=0
            for i in subset:
                #read each feature
                feat = self.layer.GetFeature(i)
                features.append(feat)

                #read each field associated to the feature
                for f in self.fields.name:
                    self.fields.values[f][k] = feat.GetField(f)
                k+=1
        self.features = np.array(features)


    def FeatureCount(self):
        """
        Return the number of features in self
        """

        return self.layer.GetFeatureCount()


    def reproject(self,target_srs):
        """
        Return a new SingleLayerVector object with features reprojected according to target_srs.
        """
        # create the CoordinateTransformation
        coordTrans = osr.CoordinateTransformation(self.srs, target_srs)

        # create the output layer
        outDataSet = ogr.GetDriverByName('Memory').CreateDataSource('Reprojected '+self.layer.GetName())
#        outputShapefile = r'c:\data\spatial\basemap_4326.shp'
#        outDataSet = driver.CreateDataSource(outputShapefile)
        outLayer = outDataSet.CreateLayer(self.layer.GetName(), target_srs,geom_type=ogr.wkbMultiPolygon)

        # add fields
        inLayerDefn = self.layer.GetLayerDefn()
        for i in range(0, inLayerDefn.GetFieldCount()):
            fieldDefn = inLayerDefn.GetFieldDefn(i)
            outLayer.CreateField(fieldDefn)

        # get the output layer's feature definition
        outLayerDefn = outLayer.GetLayerDefn()

        # loop through the input features
        self.layer.GetFeatureCount()  #to read from beginning again
        inFeature = self.layer.GetNextFeature()
        while inFeature:
            # get the input geometry
            geom = inFeature.GetGeometryRef()
            # reproject the geometry
            geom.Transform(coordTrans)
            # create a new feature
            outFeature = ogr.Feature(outLayerDefn)
            # set the geometry and attribute
            outFeature.SetGeometry(geom)
            # set the Spatial Reference
            geom.AssignSpatialReference(target_srs)
        
            for i in range(0, outLayerDefn.GetFieldCount()):
                outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
            # add the feature to the shapefile
            outLayer.CreateFeature(outFeature)
            # destroy the features and get the next input feature
            outFeature.Destroy()
            inFeature.Destroy()
            inFeature = self.layer.GetNextFeature()

        return SingleLayerVector(outDataSet)


    def draw(self,subset='all',extent='default',map_obj=None,ax='none',**kwargs):
        """
        Plot the geometries defined in the vector file
        Inputs :
        - subset : indices of the features to plot (Default is 'all')
        - extent : xmin, xmax, ymin, ymax (Default is self.extent)
        - map_obj : Basemap object, used to plot on a map
        **kwargs : any optional argument accepted by the matplotlib.patches.PathPatch class, e.g.
            - edgecolor : mpl color spec, or None for default, or ‘none’ for no color
            - facecolor : mpl color spec, or None for default, or ‘none’ for no color
            - lw : linewidth, float or None for default
            - alpha : transparency, float or None
        """

        if not hasattr(self,'features'):
            print "You must run self.read() first"
            return 0

        if subset=='all':
            subset = range(len(self.features))
        elif isinstance(subset,numbers.Number):
            subset = [subset,] #create list if only one value

        p = []
        for feat in self.features[subset]:
            sh = Shape(feat)
            if map_obj==None:
                p0 = sh.draw(ax=ax,**kwargs)
            else:
                p0 = sh.draw_on_map(map_obj,ax=ax,**kwargs)
            p.append(p0)
        
        if ax=='none':
            ax = pl.gca()
        if extent=='default':
            if map_obj==None:
                xmin, xmax,ymin,ymax = self.extent
            else:
                xmin, xmax, ymin, ymax = map_obj.xmin, map_obj.xmax, map_obj.ymin, map_obj.ymax
        else:
            xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)

        return p

    def draw_by_attr(self,attr,cmap=cm.jet,subset='all',vmin='default',vmax='default',**kwargs):
        """
        Plot the geometries defined in the vector file
        Inputs :
        - subset : indices of the features to plot (Default is 'all')
        **kwargs : any optional argument accepted by the matplotlib.patches.PathPatch class, e.g.
            - edgecolor : mpl color spec, or None for default, or ‘none’ for no color
            - facecolor : mpl color spec, or None for default, or ‘none’ for no color
            - lw : linewidth, float or None for default
            - alpha : transparency, float or None
        """

        if not hasattr(self,'features'):
            print "You must run self.read() first"
            return 0

        if subset=='all':
            subset = range(len(self.features))
        elif isinstance(subset,numbers.Number):
            subset = [subset,] #create list if only one value

        #create a collection of patches
        patches = []
        for k in subset:
            sh = Shape(self.features[k])
            patches.append(PathPatch(sh.path))

        #remove Nan values
        values = np.copy(self.fields.values[attr])
        if vmin=='default':
            vmin = np.nanmin(values)
        if vmax=='default':
            vmax = np.nanmax(values)
#        p5 = np.percentile(values[~np.isnan(values)],5)
#        p95 = np.percentile(values[~np.isnan(values)],95)
        bounds = np.linspace(vmin,vmax,255)
        norm = colors.BoundaryNorm(bounds, cmap.N)
        values[np.isnan(values)] = -1e5
        cmap = deepcopy(cmap)   #copy to avoid changes applying to other scripts
        cmap.set_under('grey')

        p = PatchCollection(patches, cmap=cmap,norm=norm)
        p.set_array(np.array(values)) #set colors

        #Plot
        ax = pl.gca()
        ax.add_collection(p)
        cb = pl.colorbar(p)
        cb.set_label(attr)
        xmin, xmax,ymin,ymax = self.extent
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)


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
        self.layer.SetSpatialFilter(poly)
        
        #update attributes
        self.extent = poly.GetEnvelope()
    

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


    def zonal_statistics(self,rasterfile,subset='all',nodata=None):
        """
        Compute statistics of the data in rasterfile for each feature in self.
        Inputs :
        - subset : indices of the features to compute (Default is 'all')
        """

        # Read raster coordinates
        img = raster.SingleBandRaster(rasterfile)
        X, Y = img.coordinates()

        # Coordinate transformation from vector to raster projection system
        coordTrans = osr.CoordinateTransformation(self.srs,img.srs)

        #Select only subset specified by user
        if subset=='all':
            subset = range(len(self.features))
        elif isinstance(subset,numbers.Number):
            subset = [subset,] #create list if only one value

        print "Loop on all features"
        median = []
        std = []
        count = []
        frac = []
        for feat in self.features[subset]:

            #Read feature geometry and reproject to raster projection
            sh = Shape(feat,load_data=False)
            sh.geom.Transform(coordTrans)
            sh.read()

            #Look only for points within the box of the feature
            xmin, xmax, ymin, ymax = sh.geom.GetEnvelope()
            inds = np.where((X>=xmin) & (X<=xmax) & (Y>=ymin) & (Y<=ymax))

            #Select data within the feature
            inside = geo.points_inside_polygon(X[inds],Y[inds],sh.vertices,skip_holes=False)
            inside = np.where(inside)
            inside_i = inds[0][inside]
            inside_j = inds[1][inside]

            data = img.r[inside_i,inside_j]

            #Filter no data values
            data = data[~np.isnan(data)]
            if nodata!=None:
                data = data[data!=nodata]

            median.append(data)

            #Compute statistics
            if len(data)>0:
                median.append(np.median(data))
                mad = 1.4826*np.median(np.abs(data-np.median(data)))
                std.append(mad)
                #                std.append(np.std(data))
                count.append(len(data))
                frac.append(float(len(data))/len(inside_i))
            else:
                median.append(np.nan)
                std.append(np.nan)
                count.append(0)
                frac.append(0)

        return np.array(median), np.array(std), np.array(count), np.array(frac)

    def create_mask(self,raster):
        """
        Return a mask (Byte) of the polygons in self, in the Spatial Reference and resolution of raster
        """

        # Open the raster file
        xsize, ysize = raster.r.shape
        x_min, x_max, y_min, y_max = raster.extent
    

        # Create memory target raster
        target_ds = gdal.GetDriverByName('MEM').Create('', ysize, xsize, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform((
                x_min, raster.xres, 0,
                y_max, 0, raster.yres,
                ))
    
        # Make the target raster have the same projection as the source raster
        target_ds.SetProjection(raster.srs.ExportToWkt())

        # Rasterize
        err = gdal.RasterizeLayer(target_ds, [1], self.layer,burn_values=[255])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)

        # Read mask raster as arrays
        bandmask = target_ds.GetRasterBand(1)
        datamask = bandmask.ReadAsArray(0, 0, ysize, xsize)

        return datamask


    def create_mask_attr(self,raster,attr):
        """
        Return a raster (Float32) containing the values of attr for each polygon in self, in the Spatial Reference and resolution of raster
        """

        # Open the raster file
        xsize, ysize = raster.r.shape
        x_min, x_max, y_min, y_max = raster.extent
    

        # Create memory target raster
        target_ds = gdal.GetDriverByName('MEM').Create('', ysize, xsize, 1, gdal.GDT_Float32)
        target_ds.SetGeoTransform((
                x_min, raster.xres, 0,
                y_max, 0, raster.yres,
                ))

        # Make the target raster have the same projection as the source raster
        target_ds.SetProjection(raster.srs.ExportToWkt())

        #loop on features
        for feat in self.features:

            # Create a memory layer to rasterize from.
            input_raster = ogr.GetDriverByName('Memory').CreateDataSource('')
            layer = input_raster.CreateLayer('poly', srs=self.srs)

            # Add polygon
    #        feat = ogr.Feature(layer.GetLayerDefn())
    #        feat.SetGeometry(feature.GetGeometryRef())
            layer.CreateFeature(feat)

            # Rasterize
            err = gdal.RasterizeLayer(target_ds, [1], layer,burn_values=[feat.GetField(attr)])

            if err != 0:
                raise Exception("error rasterizing layer: %s" % err)

            # Read mask raster as arrays
            bandmask = target_ds.GetRasterBand(1)
            datamask = bandmask.ReadAsArray(0, 0, ysize, xsize)


            # import pylab as pl
            # pl.imshow(datamask)
            # pl.colorbar()
            # pl.show()

        return datamask


    def get_extent_projected(self,pyproj_obj):
        """ Return vector extent in a projected coordinate system.

        This is particularly useful for converting vector extent to the 
        coordinate system of a Basemap instance.

        Parameters:
            pyproj_obj : A pyproj instance (such as a Basemap instance) of the
                system to convert into.

        Returns:
            (left,right,bottom,top)

        """
        xll,xur,yll,yur = self.extent

        left,bottom = pyproj_obj(xll,yll)
        right,top = pyproj_obj(xur,yur)
        return (left,right,bottom,top)




class Shape():

        def __init__(self,feature,load_data=True):
            
            self.feature = feature
            self.geom = deepcopy(feature.GetGeometryRef()) #deepcopy otherwise the feature can be modified

            if load_data==True:
                self.read()


        def draw(self, ax='none', **kwargs):
            """
            Draw the shape using the matplotlib.path.Path and matplotlib.patches.PathPatch classes
            **kwargs : any optional argument accepted by the matplotlib.patches.PathPatch class, e.g.
            - edgecolor : mpl color spec, or None for default, or ‘none’ for no color
            - facecolor : mpl color spec, or None for default, or ‘none’ for no color
            - lw : linewidth, float or None for default
            - alpha : transparency, float or None
            """
            
            #Create patch and get extent
            patch = PathPatch(self.path,**kwargs)
            xmin,xmax,ymin,ymax = self.extent

#            fig = pl.figure()
#            ax=fig.add_subplot(111)
            if ax=='none':
                ax = pl.gca()
            ax.add_patch(patch)
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)

            return patch

        def read(self):
            """
            Read shape extent and vertices
            """

            self.extent = self.geom.GetEnvelope()

            vertices = []
            codes = []

            #POLYGONS
            if self.geom.GetGeometryName()=='POLYGON':
                for i in xrange(self.geom.GetGeometryCount()):
                    poly = self.geom.GetGeometryRef(i)
                    for j in xrange(poly.GetPointCount()):
                        vertices.append([poly.GetX(j),poly.GetY(j)])
                        if j==0:
                            codes.append(Path.MOVETO)
                        elif j==poly.GetPointCount()-1:
                            codes.append(Path.CLOSEPOLY)
                        else:
                            codes.append(Path.LINETO)
        
            elif self.geom.GetGeometryName()=='MULTIPOLYGON':
                for i in xrange(self.geom.GetGeometryCount()):
                    poly = self.geom.GetGeometryRef(i)

                    for j in xrange(poly.GetGeometryCount()):
                        ring = poly.GetGeometryRef(j)
                        for k in xrange(ring.GetPointCount()):
                            vertices.append([ring.GetX(k),ring.GetY(k)])
                            if k==0:
                                codes.append(Path.MOVETO)
                            elif k==ring.GetPointCount()-1:
                                codes.append(Path.CLOSEPOLY)
                            else:
                                codes.append(Path.LINETO)
    
            #LINESTRING
            elif self.geom.GetGeometryName()=='LINESTRING':
                for j in xrange(self.geom.GetPointCount()):
                    vertices.append([self.geom.GetX(j),self.geom.GetY(j)])
                    if j==0:
                        codes.append(Path.MOVETO)
                    else:
                        codes.append(Path.LINETO)

            #MULTILINESTRING
            elif self.geom.GetGeometryName()=='MULTILINESTRING':
                for i in xrange(self.geom.GetGeometryCount()):
                    poly = self.geom.GetGeometryRef(i)
                    for k in xrange(poly.GetPointCount()):
                        vertices.append([poly.GetX(k),poly.GetY(k)])
                        if k==0:
                            codes.append(Path.MOVETO)
                        else:
                            codes.append(Path.LINETO)

            else:
                print "Geometry type %s not implemented" %self.geom.GetGeometryName()
                sys.exit(1)

            self.vertices = vertices   #list of vertices
            self.codes = codes
            self.path = Path(vertices,codes)  #used for plotting

        def draw_on_map(self,m,ax='none', **kwargs):
            """
            Draw the shape using the matplotlib.path.Path and matplotlib.patches.PathPatch classes
            **kwargs : any optional argument accepted by the matplotlib.patches.PathPatch class, e.g.
            - edgecolor : mpl color spec, or None for default, or ‘none’ for no color
            - facecolor : mpl color spec, or None for default, or ‘none’ for no color
            - lw : linewidth, float or None for default
            - alpha : transparency, float or None
            """
            
            #Convert coordinates to Basemap coordinates
            vertices = []
            for vertice in self.vertices:
                x,y = vertice
                vertices.append(m(x,y))

            #create patch
            path = Path(vertices,self.codes)
            patch = PathPatch(path,**kwargs)

            #Get extent
            xmin,xmax,ymin,ymax = self.extent
            xmin, ymin = m(xmin,ymin)
            xmax, ymax = m(xmax,ymax)

            #plot
            if ax=='none':
                ax = pl.gca()
            ax.add_patch(patch)
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)

            return patch

        def rasterize(self,srs,pixel_size,extent='default'):
            
            # Create a memory raster to rasterize into.
            if extent=='default':
                xmin, xmax, ymin, ymax = self.extent
            else:
                xmin, xmax, ymin, ymax = extent
            xsize = int((xmax-xmin)/pixel_size)
            ysize = int((ymax-ymin)/pixel_size)
            target_ds = gdal.GetDriverByName('MEM').Create( '', xsize, ysize, bands=1,eType=gdal.GDT_Byte)
            target_ds.SetGeoTransform((xmin,pixel_size,0,ymax,0,-pixel_size))
            target_ds.SetProjection(srs.ExportToWkt())

            # Create a memory layer to rasterize from.
            input_raster = ogr.GetDriverByName('Memory').CreateDataSource('')
            layer = input_raster.CreateLayer('poly', srs=srs)
            
            # Add polygon
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetGeometry(self.geom)
            layer.CreateFeature(feat)
        
            # Rasterize
            err = gdal.RasterizeLayer(target_ds, [1], layer,burn_values=[255])

            if err != 0:
                raise Exception("error rasterizing layer: %s" % err)

            # Read mask raster as arrays
            bandmask = target_ds.GetRasterBand(1)
            mask = bandmask.ReadAsArray(0, 0, xsize, ysize)

            trans = target_ds.GetGeoTransform()
            shape = (xsize, ysize)
            x = np.array(np.linspace(trans[0],(trans[0]+(xsize*trans[1])),xsize).tolist() * ysize).reshape(shape)
            y = np.array(np.linspace(trans[3],(trans[3]+(ysize*trans[5])),ysize).tolist() * xsize).reshape(shape[::-1]).T
            
            return mask, x, y


        def centerline(self,srs,pixel_size):
            """
            Create a raster of the shape and display to help user drawing the centerline
            srs : OSR SpatialReference, used for rasterization
            pixel_size : size of the raster pixel
            Outputs :
                x and y coordinates of the user clicks
            """

            # Create a binary raster of the shape
            mask, xx, yy = self.rasterize(srs,pixel_size)

            # Fill small holes on mask, e.g supraglacial lakes
            mask_fill = ndimage.binary_opening(mask,iterations=2)
            
            # Compute the skeleton of the mask to help user
            skel = morphology.skeletonize(mask_fill > 0).astype('float32')
            skel[skel==0] = np.nan

            # Show mask and skeleton
            pl.imshow(mask_fill,cmap=pl.get_cmap('Greys_r'))
            pl.imshow(skel)

            #User click to create profile
            print("Select your profile")
            coords=pl.ginput(n=0,timeout=0)
            
            #Convert to list of indices
            coords = np.int32(coords)
            coords = np.transpose(coords)
            coords = (coords[0],coords[1])

            return xx[coords], yy[coords]

                      
                      
