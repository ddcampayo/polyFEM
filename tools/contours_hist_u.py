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

import glob

#dirs=glob.glob('[0-9]*')
dirs=glob.glob('[5-9]??')

dirs.sort( key = float )

import pylab as pl

LL=1

for dir_step in dirs[1:] :

    step = int(dir_step)

 #   dtm=pl.loadtxt(path+str(n)+'/particles.dat')
    dtm=np.loadtxt(path+str(step)+'/mesh.dat')
    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]
#    rr=max(vxm)/10
#    p1 = dtm[:,5] - pl.mean(dtm[:,5] )
#    dt2=pl.loadtxt(path+'fem1/'+str(n)+'/particles.dat')
#    dt2=pl.loadtxt(path+str(n)+'/mesh.dat')

    X, Y, Z = grid(xm, ym, vxm)

#    pl.contourf( X,Y,Z, [-rr,rr] , colors='k')
#    plt.contourf( X,Y,Z, [-rr,rr] )

#    fig = plt.figure()
#    ax1  = fig.add_subplot(1,2,1, adjustable='box', aspect = 1)
#    ax2  = fig.add_subplot(1,2,2)

    fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})

    ax1.set_aspect('equal', adjustable='box')

#    ax1.quiver( xm , ym , vxm , vym )

    ax1.set_xlim([ -LL/2 , LL/2 ])
    ax1.set_ylim([ -LL/2 , LL/2 ])

    ax1.contourf( X,Y,Z )
#    plt.contour( X,Y,Z ,  3, colors='k')


    dtx = np.loadtxt(str(step)+'/histo_vel_x.dat', dtype = np.float64)

    qq=dtx[ 1: , 0]

    vx=dtx[ 1: , 1]

    dty = np.loadtxt(str(step)+'/histo_vel_y.dat', dtype = np.float64)

    vy=dty[ 1: , 1]

    # Kinetic energy
    #vv = (vx**2 + vy**2)/2

    vv = (vx + vy)/2

    tiny = 0 # 1e-16

    #ax2.semilogy( qq , vv + tiny )
    ax2.loglog( qq , vv + tiny )

    ax2.set_ylim( 1e-15 , 1e-1 )

    fig.tight_layout()

    plt.savefig('./contours_histo/contours_histo_'+str(step).zfill(5)+'.png')

    plt.close( fig )
#    pl.colorbar(ticks=[0.45,0.55])

#pl.xlim([-0.7,0.7]) ; pl.xlabel('$x/L$')
#pl.ylim([-0.7,0.7]) ; pl.ylabel('$y/L$')
