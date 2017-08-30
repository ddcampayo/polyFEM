#!/usr/bin/python

from scipy import integrate

import numpy as np
import sys

#LL=128
#NN=512

if len(sys.argv) == 1:
   exp= 1.0
else:
   exp=float(sys.argv[1])

i_exp = 1.0/exp

import pylab as pl

#pl.figure(figsize=(8,8))

Dt=0.5
skip=10
begin=skip
end=100000
limits=1
#path='timings_full/'
path='./'

#dq = 2 * pl.pi /  LL

for n in range(begin,end+skip,skip):
    pl.clf()
#    dt=pl.loadtxt(str(n)+'/histo_phi_'+str(n)+'.dat')
    dt=pl.loadtxt(str(n)+'/histo_phi.dat')

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

    L1 = 2 * pl.pi /  q1

#    print n,q1,L1
    print n*Dt , L1

#    pl.plot( n,q1)

#savefig('q1_'+str(n/skip)+'.png')
