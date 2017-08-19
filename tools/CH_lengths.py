#!/usr/bin/python
from scipy import integrate

LL=128
NN=512

import pylab as pl

pl.figure(figsize=(8,8))

skip=100
begin=0000
end=100000
limits=1
#path='timings_full/'
path='./'

#dq = 2 * pl.pi /  LL

for n in range(begin,end+skip,skip):
    pl.clf()
    dt=pl.loadtxt('histo_phi_'+str(n)+'.dat')

    qq=dt[:,0]
    ff=dt[:,1]

    #q1 = dq*sum(qq*ff)
    ############  integrate.trapz
    intS =  integrate.trapz(ff, qq )
    intqS =  integrate.trapz(qq*ff, qq )
    #    q1 = sum(qq*ff) / sum(ff)
    q1 = intqS / intS

    L1 = 2 * pl.pi /  q1

    print n,q1,L1

    pl.plot( n,q1)

savefig('q1_'+str(n/skip)+'.png')
