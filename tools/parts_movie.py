#!/usr/bin/python

LL=40

import pylab as pl

pl.figure(figsize=(8,8))

skip=10
begin=1400
end=20000
limits=1
#path='timings_full/'
path='./'

for n in range(begin,end+skip,skip):
    pl.clf()
    dtm=pl.loadtxt(path+str(n)+'/particles.dat')
#    dtm=pl.loadtxt(path+str(n)+'/mesh.dat')
    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

#    pl.scatter( xm , ym , c=alm, s=80, vmin= -limits , vmax= limits)
    pl.scatter( xm , ym , c=alm, s=80 )

    pl.xlim([ -LL/2 , LL/2 ])
    pl.ylim([ -LL/2 , LL/2 ])
#    pl.colorbar()
#    pl.colorbar(ticks=[0.45,0.55])
    pl.savefig('parts'+str(n/skip))
