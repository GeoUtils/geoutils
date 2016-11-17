import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import os

dpath = os.path.dirname(__file__)
class BiInterpolator:
    '''Bilinear interpolation in 2D. The code is modified from mpl_toolkits.basemap.interp and scipy.interpolate'''
    def __init__(self, xin, yin, datain):
        '''Setting up the interpolator.

        .. Args:

            * xin     -> Monotonic array of x coordinates
            * yin     -> Monotonic array of y coordinates
            * datain  -> 2D array corresponding to (y,x) '''

        if xin.shape[0] != datain.shape[1]:
            raise ValueError('Shapes of datain and x do not match')

        if yin.shape[0] != datain.shape[0]:
            raise ValueError('Shapes of datain and y do not match')

        if xin[-1] < xin[0]:
            raise ValueError('Array x not sorted')

        if yin[-1] < yin[0]:
            raise ValueError('Array y not sorted')

        self.xin = xin.copy()
        self.yin = yin.copy()

        delx = xin[1:] - xin[0:-1]
        dely = yin[1:] - yin[0:-1]

        if max(delx)-min(delx) < 1.e-4 and max(dely)-min(dely) < 1.e-4:
            self.regular = True
        else:
            self.regular = False

        self.xinlist = self.xin.tolist()
        self.yinlist = self.yin.tolist()
        self.nx = len(self.xinlist)
        self.ny = len(self.yinlist)
        self.zin = datain


    def __call__(self,xi,yi):
        '''Function call to actually interpolate.'''
        if xi.shape != yi.shape:
            raise ValueError('xi and yi must have same shape.')

        if self.regular:
            xcoords = (self.nx-1)*(xi-self.xin[0])/(self.xin[-1]-self.xin[0])
            ycoords = (self.ny-1)*(yi-self.yin[0])/(self.yin[-1]-self.yin[0])
        else:
            xiflat = xi.flatten()
            yiflat = yi.flatten()
            ix = (np.searchsorted(self.xin,xiflat)-1).tolist()
            iy = (np.searchsorted(self.yin,yiflat)-1).tolist()
            xiflat = xiflat.tolist()
            yiflat = yiflat.tolist()

            xin = self.xinlist
            yin = self.yinlist
                
            xcoords = []
            ycoords = []
            for n,i in enumerate(ix):
                if i < 0:
                    xcoords.append(-1)
                elif i >= self.nx-1:
                    xcoords.append(self.nx)
                else:
                    xcoords.append(float(i)+(xiflat[n]-xin[i])/(xin[i+1]-xin[i]))
            for m,j in enumerate(iy):
                if j < 0:
                    ycoords.append(-1)
                elif j >= self.ny-1:
                    ycoords.append(self.ny)
                else:
                    ycoords.append(float(j)+(yiflat[m]-yin[j])/(yin[j+1]-yin[j]))

            xcoords = np.reshape(xcoords, xi.shape)
            ycoords = np.reshape(ycoords, yi.shape)

        xcoords = np.clip(xcoords,0,self.nx-1)
        ycoords = np.clip(ycoords,0,self.ny-1)

        xint = xcoords.astype(np.int32)
        yint = ycoords.astype(np.int32)
        xip1 = np.clip(xint+1,0,self.nx-1)
        yip1 = np.clip(yint+1,0,self.ny-1)

        delx = xcoords - xint.astype(np.float32)
        dely = ycoords - yint.astype(np.float32)

        zin = self.zin
        dataout = (1.-delx)*(1.-dely)*zin[yint,xint] + \
                    delx*dely*zin[yip1,xip1] + \
                    (1.-delx)*dely*zin[yip1,xint] + \
                    delx*(1.-dely)*zin[yint,xip1]

        return dataout




class EGM96reader:
    '''Create a default reader. Can use it with lon/ lat as inputs.'''
    def __init__(self):
        #####Read the metadata for the raster file
        egm96full = gdal.Open('%s/egm96-5.pgm'%(dpath))
        mdict = egm96full.GetMetadata_Dict()
        zscale = np.float(mdict['Scale'])
        zoff = np.float(mdict['Offset'])

        #########Get the geo-coordinates
        (x0,delx,dum0,y0,dum1,dely)  = egm96full.GetGeoTransform()
        xsize = egm96full.RasterXSize
        ysize = egm96full.RasterYSize


        ########Note that we set up data for (lon, -lat)
        xin = x0 + delx*np.arange(xsize)
        yin = y0 + dely*np.arange(ysize)
        yin = -1.0*yin                     #Order reversed

        zin = egm96full.ReadAsArray(0,0,xsize,ysize)
        zin =zin*zscale + zoff
 
        self.oper  = BiInterpolator(xin, yin, zin)
        

    def __call__(self, xi, yi):
        #Call the object with lon and lat
        yinp = -yi
        xinp = xi.copy()
        xinp[xinp<0] += 360.0

        return self.oper(xinp,yinp)


    def test(self):
        Datin = np.loadtxt('%s/GeoidHeights.dat'%(dpath))
        lat = -1.0*Datin[:,0]    #Multiply by -1 
        lon = Datin[:,1]
        lon[lon<0] +=360.0
        testh = Datin[:,3]   

        res = self.oper(lon,lat)
        err = res-testh
        print 'All errors are reported in meters.'
        print '**********************************'
        print 'Maximum Error: ', np.max(err)
        print 'Minimum Error: ', np.min(err)
        print 'Standard deviation: ', np.std(err)


    def plot(self):
        plt.imshow(self.oper.zin,extent=[0,360,-90,90])
        plt.colorbar()
        plt.title('Geoid Offset in meters')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.show()

#   res = reader(360.0+np.array([-118.0,-118.0]), -1.0*np.array([35.0,35.0]))
#   print res
