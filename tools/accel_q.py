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

dtx = np.loadtxt(step+'/histo_accel_x.dat', dtype = np.float64)

qq=dtx[ 1: , 0]

ax=dtx[ 1: , 1]

dty = np.loadtxt(step+'/histo_accel_y.dat', dtype = np.float64)

ay=dty[ 1: , 1]

aa = ax**2 + ay**2

tiny=1e-30

pl.plot( qq , np.log( aa + tiny) )

pl.savefig('accel_'+step+'.png')
