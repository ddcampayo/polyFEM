#!/usr/bin/python
import glob

import numpy as np

from numpy import pi

Dt=0.01
nu=0.001
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

	time=Dt*step # + Dt/1

	A = (1-np.exp( - 4 * pi**2 * nu * time))/nu

	dt = np.loadtxt(dir_step+'/particles.dat')
	#dt = np.loadtxt(dir_step+'/mesh.dat')
#        x=dt[:,18]
#        y=dt[:,19]

	x=dt[:,0]
	y=dt[:,1]

	vx=dt[:,8]
	vy=dt[:,9]

        vol=dt[:,3]

	ux=  A * np.cos( 2*pi*y / LL)

	ddx=(vx-ux)**2

	dd= vol * ddx

	ff= vol * ux**2

        eps=1e-12

	print " %g  %g " % (
	time ,
	( ( eps + np.sum( dd ) ) / (eps + np.sum( ff ) ) )**0.5
	)

