.. _get-started:

Getting Started
---------------


Vector files handling within Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
   
   from geoutils import geovector as vect

Read a vector \file::
  
  outlines = vect.SingleLayerVector('data/glacier_outlines.shp')

At this stage, only the attributes are read, not the individual features. Metadata are saved as attributes, e.g.:

- extent::
    
    outlines.extent
    
- osr.SpatialReference::
    
    outlines.srs
    
- number of features selected::
    
    outlines.FeatureCount()

Actually read the features:::
  
  outlines.read()

The features are then stored as a list in outlines.features.

Display the layer::
  
  import matplotlib.pyplot as plt
  outlines.draw()
  plt.show()

Display only the contours::
  
  outlines.draw(facecolor='none')
  plt.show()

Display only a subset, e.g. glaciers larger than 5 km2::
  
  import numpy as np
  subset = np.where(outlines.fields.values['Area']>5)
  outlines.draw(facecolor='none',subset=subset)
  plt.show()

A specific outline::
  
  outlines.draw(facecolor='none',subset=np.where(outlines.fields.values['RGIId']==b'RGI50-15.09991'))
  plt.show()

Crop the layer
==============

With specified coordinates::
  
  outlines.crop(86.75,86.93,27.9,28.05)
  outlines.read()

From an external raster::

  outlines.crop2raster('data/LE71400412000304SGS00_B4_crop.TIF')
  outlines.read()


Reproject the layer
===================

::

   import georaster as raster
   img = raster.SingleBandRaster('data/LE71400412000304SGS00_B4_crop.TIF')
   outlines_proj = outlines.reproject(img.srs)
   outlines_proj.read()

Display on top of the raster::
  
  plt.imshow(img.r,cmap='gray',extent=img.extent)
  outlines_proj.draw(facecolor='none',edgecolor='r')
  plt.show()

Save the layer to an ESRI shapefile
===================================

For some reasons, the RGIId attribute does not save properly, so needs to be removed first::

  del outlines_proj.fields.values['RGIId']
  
Save to \file::

   vect.save_shapefile('temp.shp',outlines_proj)


Create a mask
=============

Mask inside/outside::
  
  ice_mask = outlines.create_mask(img)  # ice_mask=0 outside, 255 inside
  plt.imshow(ice_mask)
  plt.colorbar()
  plt.show()
  
Mask by feature attributes::
  
  area = outlines.create_mask_attr(img,'Area')
  plt.imshow(area)
  plt.colorbar()
  plt.show()

  
Clip a raster to a vector layer
===============================

Cropping but no masking outside the mask::
  
  feat=11
  clipped_data = outlines.clip_raster('data/LE71400412000304SGS00_B4_crop.TIF','none',feature=feat)
  clipped_data = raster.SingleBandRaster(clipped_data)
  plt.imshow(clipped_data.r,cmap='gray',extent=clipped_data.extent)
  outlines_proj.draw(facecolor='none',edgecolor='r',subset=feat)
  plt.show()

With masking (bugs, not only selected feature displayed + no data = 241...)::
  
  clipped_data = outlines.clip_raster('data/LE71400412000304SGS00_B4_crop.TIF','none',feature=feat,masking=True)
  clipped_data = raster.SingleBandRaster(clipped_data)
  plt.imshow(clipped_data.r,cmap='gray',extent=clipped_data.extent)
  outlines_proj.draw(facecolor='none',edgecolor='r',subset=feat)
  plt.show()

  
Zonal statistics
================

Calculate the mean of the raster inside the features::
  
  mean = outlines.zonal_statistics(img,np.mean)

Save as one of the field attributes::

  outlines.fields.values['mean'] = mean

Display the results::
  outlines.draw_by_attr('mean',cmap=plt.get_cmap('gray'))
  plt.show()


Extract raster values along profiles
====================================

Load the centerline dataset, reproject and crop to raster extent::

   cfl = vect.SingleLayerVector('data/centerlines.shp')
   cfl = cfl.reproject(img.srs)
   cfl.crop2raster(img.ds_file)
   cfl.read()

Display both together::

  plt.imshow(img.r,cmap='gray',extent=img.extent)
  cfl.draw(facecolor='none',edgecolor='r')
  plt.show()
  
Extract the raster values at specific centerlines::

  cfl.layer.SetAttributeFilter("FID=9 OR FID=73")
  xx, yy, profiles = cfl.extract_value_from_raster(img)

Plot the outputs::
  
  for p in profiles:
    plt.plot(p)
  plt.xlabel('Pixel along line')
  plt.ylabel('Image intensity')
  plt.show()

  
Command line tools
~~~~~~~~~~~~~~~~~~

Display a raster image
======================

::
   
    plot_geo_image.py data/LE71400412000304SGS00_B4_crop.TIF

Choose a different color scale, remove colorbar::
  
    plot_geo_image.py LE71400412000304SGS00_B4_crop.TIF -cmap gray -nocb

Add some labelling::
  
  plot_geo_image.py LE71400412000304SGS00_B4_crop.TIF -clabel 'Image intensity' -title 'Landsat image of the Everest region' -cmap gray

Change figure size and save to \file::
  
  plot_geo_image.py LE71400412000304SGS00_B4_crop.TIF -cmap gray -figsize 12,8 -dpi 150 -save temp.png
