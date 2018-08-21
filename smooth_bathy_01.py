from __future__ import division
#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2017, NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"


"""
Try to look for potential un-realistic wave forces

"""

import matplotlib as mpl
mpl.use('Agg')

from pyproj import Proj

import os,sys
import netCDF4 as n4
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.tri as Tri
import scipy as sp

sys.path.append('/disks/NASARCHIVE/saeed_moghimi/opt/python-packages/')
import pynmd.models.adcirc.post as adcp

#from dateutil import parser
#import pynmd.models.adcirc.post as adcp
import pynmd.plotting.colormaps as cmaps


from matplotlib.tri import Triangulation, LinearTriInterpolator
import math


def find_neighbors(pindex, triang):
    neighbors = list()
    for simplex in triang.vertices:
        if pindex in simplex:
            neighbors.extend([simplex[i] for i in range(len(simplex)) if simplex[i] != pindex])
            '''
            this is a one liner for if a simplex contains the point we`re interested in,
            extend the neighbors list by appending all the *other* point indices in the simplex
            '''
    #now we just have to strip out all the dulicate indices and return the neighbors list:
    return list(set(neighbors))

def find_nearest1d(xvec,yvec,xp,yp,n):
    """
    In: xvec, yvec of the grid and xp,yp of the point of interst
    Retun: i,j,proximity of the nearset grid point
    np : numper of closeset points
    """
    dist = np.sqrt((xvec-xp)**2+(yvec-yp)**2)
    ind = np.argsort(dist)[1:n+1]
    return ind,dist[ind]


def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude
    Not right for some polar regions (Norway, Svalbard, Antartica)
    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return np.floor( ( lon + 180 ) / 6) + 1
    
    

vv = False

fname = 'maxele.63.nc'

ncf = n4.Dataset(fname)
ncv_adc = ncf.variables
x = ncv_adc['x'][:]
y = ncv_adc['y'][:]
dep = ncv_adc['depth'][:]
el  = ncv_adc['element'][:] - 1
ncf.close()

utm_txt = str (int((utm_from_lon(((x.min()+x.max())/2)))))
myProj = Proj("+proj=utm +zone="+utm_txt+"K, +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
UTMx, UTMy = myProj(x, y)

print ' > Gen tri .. '

tri     = Tri.Triangulation(UTMx,UTMy, triangles=el)
trip    = Tri.Triangulation(x,y      , triangles=el)

sp_neibor = False
if sp_neibor:
    print ' > Gen tri_sp .. '
    tri_sp  = sp.spatial.Delaunay(np.array([[x1,y1] for x1,y1 in zip(x, y)]))


# Mask off unwanted triangles.
dep_mid = dep[tri.triangles].mean(axis=1)
mask = np.where(dep_mid < -10, 1, 0)
tri.set_mask(mask)
#####################
nitr = 5

for itr in range(nitr):
    nfname = 'dep_smooth_itr_' + str(100+itr) + '.nc'
    os.system('cp -f ' + fname + '  ' + nfname )


os.system('echo  Log smoothing &> log.txt  ')

for itr in range(nitr):
    print '   > Calculate gradient ..',itr,nitr
    os.system('echo  ' + str(itr) + ' from '+ str(nitr) + ' >> log.txt  ')
    os.system('echo  "************************************" >> log.txt  ')

    tci = LinearTriInterpolator(tri,dep)
    (ddep_dx,ddep_dy) = tci.gradient(UTMx,UTMy)
    grad_dep = np.sqrt(ddep_dx**2 + ddep_dy**2)
    
    grad_lim = 0.35
    [ind] = np.where( grad_dep > grad_dep.max() * grad_lim)
    
    #######
    if vv:
        fig, ax = adcp.make_map(res = 'm')
        vmin = 0 
        vmax = 0.5
        dv = (vmax-vmin)/50.0
        levels=np.arange(vmin,vmax+dv,dv)
        cf1  = ax.tricontourf (trip, grad_dep ,levels=levels,cmap = cmaps.jetMinWi, extend='both')#,alpha=0.8)#extend='max' )  #,extend='both'
        plt.colorbar(cf1,shrink = 0.5,ticks = [vmin,(vmin+vmax)/2,vmax])      
        ax.set_xlim(-99,-52.8)
        ax.set_ylim(5,46.3)
        ax.plot(x[ind],y[ind],'ro',ms=3)    
        plt.savefig(str(100+itr)+'_itr_dep_grad.png',dpi=450)
    
    dep_smt = np.zeros_like(dep) + dep
    
    print '   > Find nodes ..'

    str_log = ' >>>>> ITR > '   + str(int(itr)) +'  out of  ' + str(int(nitr))  +  '   >>>>> Total nodes to be smoothed > '  +str(len(ind))
    os.system('echo  ' + str_log +' >> log.txt  ')
    print str_log 
    all_nodes = {}
    for i in range(len(ind)):
        if dep_smt[ind[i]] < 250:
            if sp_neibor:
                nei_list = find_neighbors(pindex = ind[i] , triang = tri_sp)
                xn = x[nei_list]
                yn = y[nei_list]
                dist = np.sqrt ( (xn-x[ind[i]])**2 + (yn-y[ind[i]])**2) 
        
            else:
                nei_list,dist =  find_nearest1d(xvec = x , yvec = y , xp = x[ind[i]] , yp = y[ind[i]] , n = 15 )
            
            weights  = (1/(dist)) / (1/dist).sum()
            
            #wieghted
            mean_dep = (dep_smt[nei_list] * weights).sum()
            
            #simple avg
            #mean_dep = np.sum(dep_smt[nei_list])/len(nei_list)
            str_log =  str(i)+ '  ' + str(len(ind))+ '  ' + str(dep_smt[ind[i]]) + '  ' + str(mean_dep)
            print '     > ' + str_log
            os.system('echo  ' + str_log +' >> log.txt  ')
            
            all_nodes[str(i)] = dict(nei_list= nei_list,  weights = weights,dep = dep_smt,mean_dep = mean_dep  )
            
            dep_smt[ind[i]] = mean_dep
    


    nfname = 'dep_smooth_itr_' + str(100+itr) + '.nc'
    nc0    = n4.Dataset(nfname,'r+')
    ncv0   = nc0.variables
    ncv0['depth'][:] = dep_smt
    nc0.close()

    import pickle
    filehandler = open(b"dep_smt.pickle","wb")
    pickle.dump(all_nodes ,filehandler)
     
    dep = dep_smt
    del dep_smt
    
    
