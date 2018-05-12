#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import griddata

def grid(x, y, z , resX=90, resY=90):
    "Convert 3 column data to matplotlib grid"

    X, Y = np.mgrid[ min(x):max(x):(1j * resX) , min(y):max(y):(1j *resY) ]
    Z = griddata((x,y), z, (X, Y), method='linear')

    return X, Y, Z

plt.figure(figsize=(8,8))

skip=1
#path='timings_full/'
path='./'

for n in range(1,400+skip,skip):
 #   dtm=pl.loadtxt(path+str(n)+'/particles.dat')
    dtm=np.loadtxt(path+str(n)+'/mesh.dat')
    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]
#    rr=max(vxm)/10
#    p1 = dtm[:,5] - pl.mean(dtm[:,5] )
#    dt2=pl.loadtxt(path+'fem1/'+str(n)+'/particles.dat')
#    dt2=pl.loadtxt(path+str(n)+'/mesh.dat')

    X, Y, Z = grid(xm, ym, vxm)

#    pl.contourf( X,Y,Z, [-rr,rr] , colors='k')
#    plt.contourf( X,Y,Z, [-rr,rr] )

    plt.contourf( X,Y,Z )
#    plt.contour( X,Y,Z ,  3, colors='k')
    
    plt.savefig( 'contours'+str(n)+'.png' )

#    pl.colorbar(ticks=[0.45,0.55])

#pl.xlim([-0.7,0.7]) ; pl.xlabel('$x/L$')
#pl.ylim([-0.7,0.7]) ; pl.ylabel('$y/L$')
