.. _get-started:

Getting Started
---------------


Display a raster image
~~~~~~~~~~~~~~~~~~~~~~

::
   plot_geo_image.py data/LE71400412000304SGS00_B4_crop.TIF

Choose a different colorbar, remove colorbar::
  plot_geo_image.py LE71400412000304SGS00_B4_crop.TIF -cmap gray -nocb

Add some labelling::
  plot_geo_image.py LE71400412000304SGS00_B4_crop.TIF -clabel 'Image intensity' -title 'Landsat image of the Everest region' -cmap gray

Change figure size and save to \file::
  plot_geo_image.py LE71400412000304SGS00_B4_crop.TIF -cmap gray -figsize 12,8 -dpi 150 -save temp.png


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
  outlines.draw(facecolor='none',subset=np.where(outlines.fields.values['RGIId']=='RGI50-15.09991'))
  plt.show()


