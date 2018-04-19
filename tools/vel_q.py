#!/usr/bin/python

from scipy import integrate

import numpy as np
import sys

from numpy import pi

#LL=32

Dt=0.01

for line in open("simu.cfg"):
 if "step" in line:
  Dt=float(line.split()[1])

if len(sys.argv) == 1:
   step='0'
else:
   step=sys.argv[1]

import pylab as pl

dtx = np.loadtxt(step+'/histo_vel_x.dat', dtype = np.float64)

qq=dtx[ 1: , 0]

vx=dtx[ 1: , 1]

dty = np.loadtxt(step+'/histo_vel_y.dat', dtype = np.float64)

vy=dty[ 1: , 1]

# Kinetic energy
vv = (vx**2 + vy**2)/2

tiny=1e-30

pl.plot( qq , np.log( vv + tiny) )

pl.savefig('vel_'+step+'.png')

np.savetxt('power_K' + step+ '.dat', np.column_stack(( qq , vv )) )
