#!/usr/bin/python

# A test of the velocity field of a TG vortices simulation

import glob

import numpy as np

from numpy import pi

Dt=0.01
nu=0.01
LL=1

for line in open("simu.cfg"):
 if "viscosity" in line:
  nu=float(line.split()[1])
 if "step" in line:
  Dt=float(line.split()[1])

dirs=glob.glob('[0-9]*')

#dirs = sorted(dirs)

dirs.sort( key = float )

for dir_step in dirs :

        step = int(dir_step)

        time=Dt*step # + Dt/1.0
        
        A = np.exp( -8 * pi**2 * nu * time)

#        print(' t = ' , time,  ' A = ' , A )

        #dt = np.loadtxt(dir_step+'/particles.dat', dtype = np.float64)
        dt = np.loadtxt(dir_step+'/mesh.dat', dtype = np.float64)
        #        x=dt[:,18]
        #        y=dt[:,19]

        x=dt[:,0]
        y=dt[:,1]

        vx=dt[:,8]
        vy=dt[:,9]

#        vx=dt[:,15]
#        vy=dt[:,16]
        
        vol=dt[:,3]

        # vol=1

        fx=  A*np.sin( 2*pi * x / LL )*np.cos( 2*pi * y / LL)
        fy=- A*np.cos( 2*pi * x / LL )*np.sin( 2*pi * y / LL)

        ddx=(vx-fx)**2
        ddy=(vy-fy)**2
        #       ddx=np.fabs(vx-fx)
        #       ddy=np.fabs(vy-fy)

        dd= vol * (ddx+ddy)

#        print( 'dd = ' , np.sum(ddx) )
        
        ff= vol * (fx**2 + fy**2)

#       ff= vol * (np.fabs(fx) + np.fabs(fy) )


        print ( ' %g %g ' % ( time , ( np.sum( dd )  / np.sum( ff ) )**0.5 ) )

