#/usr/bin/python
# coding: utf-8

"""
Module description : Module containing useful functions for handling polygons
Author : Amaury Dehecq
"""


import numpy as np
from math import *



def unique(ar):
    """
    Returns the indices of ar that result in a unique array of pairs
    Input : ar, numpy array of (x,y) coordinates
    Output : indices, such that ar[indices] is unique
    """
    #sort pairs by y then by x
    order = np.lexsort(ar.T)
    ar = ar[order]
    #Compute the x and y difference between 2 neighbors
#    diff = ar - np.roll(ar,-1,axis=0)  #will return the last occurence
    diff = ar - np.roll(ar,1,axis=0)   #will return the 1st occurence
    #Points where any of the difference is different from 0 is unique
    inds_sorted = np.where(np.any(diff!=0,axis=1))[0]
    #convert to index of the original array, not sorted
    inds = order[inds_sorted]
    return np.sort(inds)


def point_inside_polygon(x,y,poly):
    """
    Determine if a point is inside a given polygon or not
    x,y : f, coordinates of the point to test
    poly : list or np array,  list of (x,y) pairs defining the vertex of the polygon

    Algorithm : this function counts the number of times that the half-line (D) parallel to the x-axis and reaching (x,y) from the right crosses a side of the polygon
    """
    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    #Loop on all sides of the polygon (P1,P2)
    for i in range(n+1):
        p2x,p2y = poly[i % n]

        #if the next 2 conditions not verified, (D) won't cross (P1,P2)
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):

                #If next condition is not verified, (D) coming from the right, will not cross (P1,P2)
                if x <= max(p1x,p2x):

                    if p1y != p2y:
                        #if p1y=p2y, xinters has an infinite number of solutions
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    
                    #We count intersection for the half-line coming from the right, all intersection points to the left of (x,y) are not taken into account
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside


def points_inside_polygon(x,y,poly,skip_holes=False):
    """
    Same as point_inside_polygon but x,y can be numpy.arrays containing multiple points
    Return a boolean numpy.array
    """

    n = len(poly)
    inside = np.zeros(x.shape,dtype=int)
    
    p1x,p1y = poly[0]

    if skip_holes==True:
        vertex = unique(poly)
    else:
        vertex = range(n)

    #Loop on all sides of the polygon (P1,P2)
    for i in vertex:
        p1x,p1y = poly[i % n]        
        p2x,p2y = poly[(i+1) % n]

        candidates = np.where((y > min(p1y,p2y)) & (y <= max(p1y,p2y)) & (x <= max(p1x,p2x)))[0]
        
        if len(candidates!=0):
            if p1y != p2y:
                xinters = (y[candidates]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    
                if p1x == p2x:
                    inside[candidates] = 1-inside[candidates]
                else:
                    inside[candidates] = np.where(x[candidates]<=xinters,1-inside[candidates],inside[candidates])
    
#        p1x,p1y = p2x,p2y

    inside = np.bool8(inside)
    return inside



def poly_area(lon,lat):
    """
    Compute area of a polygon on Earth
    lon, lat : np array, longitudes and latitudes of the points of the polygon
    """

    #Conversion en radian
    lon = pi*lon/180.
    lat = pi*lat/180.

    n = len(lon)
    s=0
    for i in xrange(-1,n-1):
        s+=lon[i]*lat[i+1] - lon[i+1]*lat[i]
    
    Rt = 6371
    s = s*Rt**2*cos(np.mean(lat))
    return 0.5*abs(s)


def dist_ortho(lon1,lat1,lon2,lat2):
    """
    Compute the orthodromic distance between points (lon1,lat1) and (lon2,lat2)
    Coordinates are given in decimal degrees
    One of the pair of coordinates can be an array
    """

    Re = 6378137.0    #WGS-84 equatorial radius
    
    return Re*np.arccos(np.cos(lat1*pi/180)*np.cos(lat2*pi/180)*np.cos((lon2-lon1)*pi/180) + np.sin(lat1*pi/180)*np.sin(lat2*pi/180))


class Line():

    def __init__(self,pts1,pts2):
        """
        Create a line with equation y=ax+b from 2 points on the line
        """
        
        x1, y1 = pts1
        x2, y2 = pts2

        a = (y2-y1)/(x2-x1)
        b = y1 - a*x1

        self.a = a
        self.b = b


        
    def __call__(self,x):
        """
        Return the y value on the line with absciss x
        """
        return self.a*x + self.b

    
    def within_circle(self,c,R):
        """
        Return all points of the line within the circle of center c and radius R
        """

        Xc, Yc = c
        X0 = np.arange(Xc-R,Xc+R)
        Y0 = self(X0)

        dist = np.sqrt((X0-Xc)**2+(Y0-Yc)**2)
        valid_inds = np.where(dist<=R)
        
        return X0[valid_inds], Y0[valid_inds]
        
