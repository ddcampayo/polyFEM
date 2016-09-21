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

LL=1

dirs=glob.glob('[0-9]*')

#dirs = sorted(dirs)

dirs.sort( key = float )

for dir_step in dirs :

	step = int(dir_step)

	time = Dt*step # + Dt/2.0

	A = np.exp( -16 * pi**2 * nu * time)

#	dt = np.loadtxt(dir_step+'/particles.dat')
	dt = np.loadtxt(dir_step+'/mesh.dat')
	x=dt[:,0]
	y=dt[:,1]
#	x=dt[:,18]
#	y=dt[:,19]
	p=dt[:,5]

	p0=  A/4*(np.cos(4*pi*x/LL)+np.cos(4*pi*y/LL))
        p_off= p0[0] - p[0]

#        p0 -=  p_off

        p=p-np.average(p)
        p0=p0-np.average(p0)

	dd=(1.0*p-p0)**2

	ff= p0**2

#	print " %g  %g %i" % (
	print " %g  %g " % (
	time ,
	(np.sum( dd )  / np.sum( ff ))**1
	)

