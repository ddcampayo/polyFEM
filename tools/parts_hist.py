#!/usr/bin/python

from scipy import integrate

import numpy as np
import sys

from numpy import pi

LL=1.0

Dt=0.01

for line in open("simu.cfg"):
 if "step" in line:
  Dt=float(line.split()[1])


import pylab as pl

import glob

dirs=glob.glob('[0-9]*')

dirs.sort( key = float )

for dir_step in dirs[1:] :

    step = int(dir_step)

    time=Dt*step #+ Dt/2.0
    pl.clf()

    dtm=np.loadtxt(str(step)+'/mesh.dat')
    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];
    vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

    vel=vxm**2+vym**2

    fig = pl.figure()

    ax1  = fig.add_subplot(2,1,1, adjustable='box', aspect = 1)

    ax2  = fig.add_subplot(2,1,2)

#    ax1.set_title('Sharing Y axis')

    #    pl.scatter( xm , ym , c=alm, s=80, vmin= -limits , vmax= limits)
#    ax1.scatter( xm , ym , c= pm , s=8/LL , cmap= pl.cm.bwr , linewidths=0 )

    ax1.scatter( xm , ym , c = vel , s=8/LL , cmap= pl.cm.bwr , linewidths=0 )
    
#    ax1.colorbar()
    ax1.quiver( xm , ym , vxm , vym )

    ax1.set_xlim([ -LL/2 , LL/2 ])
    ax1.set_ylim([ -LL/2 , LL/2 ])
    #    pl.colorbar(ticks=[0.45,0.55])
    #    pl.savefig('parts'+str(n/skip)  )

    dtx = np.loadtxt(str(step)+'/histo_vel_x.dat', dtype = np.float64)

    qq=dtx[ 1: , 0]

    vx=dtx[ 1: , 1]

    dty = np.loadtxt(str(step)+'/histo_vel_y.dat', dtype = np.float64)

    vy=dty[ 1: , 1]

    # Kinetic energy
    #vv = (vx**2 + vy**2)/2

    #I _think_ velocities are already squared

    vv = (vx + vy)/2

    tiny=1e-16

    ax2.semilogy( qq , vv + tiny )

    pl.savefig('parts_histo_'+str(step).zfill(5)+'.png')


    
#np.savetxt('power_K' + step+ '.dat', np.column_stack(( qq , vv )) )



