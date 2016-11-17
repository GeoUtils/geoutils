"""
Test EGM96 module
"""

import EGM96
import numpy as np

h = EGM96.EGM96reader()

h.plot()   #To plot global offsets
#h.test()   # To test on calibrated data 

#Compute offset for a region
lat = np.arange(100)*0.01+35.0
lon = np.arange(100)*0.01-118.0
offset = h(lon, lat)

print offset
