#!/usr/bin/python

from scipy import integrate

import numpy as np
import sys

from numpy import pi

LL=128

Dt=0.01

for line in open("simu.cfg"):
 if "step" in line:
  Dt=float(line.split()[1])

import glob

dirs=glob.glob('[0-9]*')

dirs.sort( key = float )


if len(sys.argv) == 1:
   exp= 1.0
else:
   exp=float(sys.argv[1])

i_exp = 1.0/exp

import pylab as pl


for dir_step in dirs[1:] :

#dq = 2 * pl.pi /  LL

#    pl.clf()
#    dt=pl.loadtxt(str(n)+'/histo_phi_'+str(n)+'.dat')
 step = int(dir_step)

 time=Dt*step #+ Dt/2.0

 dt = np.loadtxt(dir_step+'/histo_phi.dat', dtype = np.float64)

 qq=dt[ 1: , 0]
    # qq[0] = 1e-10
 ff=dt[ 1: , 1]

    #q1 = dq*sum(qq*ff)
    ############  integrate.trapz
    #intS =  integrate.trapz( ff / qq , qq )
    #intqS =  integrate.trapz( np.power( qq , exp - 1 ) * ff, qq )

 intS =  integrate.trapz( ff , qq )
 intqS =  integrate.trapz( np.power( qq , exp ) * ff, qq )

 #    q1 = sum(qq*ff) / sum(ff)

 q1 = np.power ( intqS / intS , i_exp )

 L1 = 2 * pi /  q1

#    print n,q1,L1
 print time , L1

#    pl.plot( n,q1)

#savefig('q1_'+str(n/skip)+'.png')
