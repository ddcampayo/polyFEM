#!/usr/bin/python

import pylab as pl

pl.figure(figsize=(8,8))

skip=10
#path='timings_full/'
path='./'

for n in range(0,400+skip,skip):
    pl.clf()
    dtm=pl.loadtxt(path+str(n)+'/particles.dat')
#    dtm=pl.loadtxt(path+str(n)+'/mesh.dat')
    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

    pl.scatter( xm , ym , c=alm, s=20, vmin=0, vmax=1)

    pl.xlim([-1.5,1.5])
    pl.ylim([-1.5,1.5])
#    pl.colorbar(ticks=[0.45,0.55])
    pl.savefig('parts'+str(n/skip))