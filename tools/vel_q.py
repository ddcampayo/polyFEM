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

qq = dtx[ 1: , 0]

vx=dtx[ 1: , 1]

dty = np.loadtxt(step+'/histo_vel_y.dat', dtype = np.float64)

vy=dty[ 1: , 1]

# Kinetic energy
#vv = (vx**2 + vy**2)/2

#I _think_ velocities are already squared

vv = (vx + vy)/2

Etot =  integrate.trapz( vv , qq )

# Normalize to 1
vv /= Etot

tiny=0 # 1e-16

pl.plot( qq , np.log( vv + tiny) )

pl.savefig('vel_'+step+'.png')

# print only non-vanishing values:

#mat = np.column_stack( ( qq , vv ) )

outf = open('power_K' + step+ '.dat','w')

idx = 0

for index, vv2 in np.ndenumerate(vv) :
   if (vv2 > tiny) :
     outf.write( ' %g %g \n' % ( qq[idx] , vv2 ) )

   idx += 1

outf.close()

