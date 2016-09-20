#!/usr/bin/python
import glob

import numpy as np

from numpy import pi

Dt=0.01
nu=0.001

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

	time=Dt*step # + Dt/1

	A = np.exp( -2 * pi**2 * nu * time)

	dt = np.loadtxt(dir_step+'/particles.dat')
	#dt = np.loadtxt(dir_step+'/mesh.dat')
#        x=dt[:,18]
#        y=dt[:,19]

	x=dt[:,0]
	y=dt[:,1]

	vx=dt[:,8]
	vy=dt[:,9]

	fx=  A*np.sin(pi*x)*np.cos(pi*y)
	fy=- A*np.cos(pi*x)*np.sin(pi*y)

	ddx=(vx-fx)**2
	ddy=(vy-fy)**2

	dd=ddx+ddy

	ff= fx**2 + fy**2

	print " %g  %g " % (
	time ,
	(np.sum( dd )  / np.sum( ff ))**1
	)

