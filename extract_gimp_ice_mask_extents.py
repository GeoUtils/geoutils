#!/usr/bin/env python
import os
import numpy as np

import georaster

rows = np.arange(0,6)
cols = np.arange(0,6)

L0DIR = '/disk/scratch/local.2/L0data/GIMP/'
filestart = 'GimpIceMask_15m_tile'

with open('GimpIceMask_15m_extents.csv','w') as f_out:
	f_out.write('row,col,xmin,xmax,ymin,ymax\n')
	for r in rows:
		for c in cols:
			print r,c
			f = filestart + str(r) + '_' + str(c) + '.tif'
			im = georaster.SingleBandRaster(L0DIR + f)
			ex = im.extent
			f_out.write(str(r) + ',' + str(c) + ',' + str(ex[0]) + ',' + \
				str(ex[1]) + ',' + str(ex[2]) + ',' +str(ex[3]) + '\n')
f_out.close()





