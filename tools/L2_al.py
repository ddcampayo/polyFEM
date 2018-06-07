#!/usr/bin/python
import glob

import numpy as np

from numpy import pi

Dt=0.01
nu=0.01
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

 time=Dt*(step+0)

 A = np.exp( -  2 * pi**2 * nu * time)

# dt = np.loadtxt(dir_step+'/particles.dat')
 dt = np.loadtxt(dir_step+'/mesh.dat')

# if step==0:
 x=dt[:,0]
 y=dt[:,1]
# else:
#  x=dt[:,18]
#  y=dt[:,19]


 p=dt[:,4]

 p0=  A*(np.sin(pi*x)*np.sin(pi*y))

#        p=p-np.average(p)
#        p0=p0-np.average(p0)

 dd=(p-p0)**2

 ff= p0**2

 print( " %g  %g " % (
  time ,
  (np.sum( dd )  / np.sum( ff ))**1
  )
  )

