#!/usr/bin/python

import pylab as pl

def grid(x, y, z , resX=100, resY=100):
    "Convert 3 column data to matplotlib grid"
    xi = pl.linspace(min(x), max(x), resX)
    yi = pl.linspace(min(y), max(y), resY)
    Z = pl.griddata(x, y, z, xi, yi , interp='linear')
    X, Y = pl.meshgrid(xi, yi )
    return X, Y, Z

pl.figure(figsize=(8,8))

skip=1
#path='timings_full/'
path='./'

LL=32.0

for n in range(0,2000000+skip,skip):
    pl.clf()
    dtm=pl.loadtxt(path+str(n)+'/particles.dat')
#    dtm=pl.loadtxt(path+str(n)+'/mesh.dat')
    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]
#    p1 = dtm[:,5] - pl.mean(dtm[:,5] )
#    dt2=pl.loadtxt(path+'fem1/'+str(n)+'/particles.dat')
#    dt2=pl.loadtxt(path+str(n)+'/mesh.dat')
    X, Y, Z = grid(xm, ym, alm)

    pl.contourf( X,Y,Z,[-0.2,-0.1,0,0.1,0.2])

    pl.xlim([-LL/2.0 , LL/2.0 ])
    pl.ylim([-LL/2.0 , LL/2.0 ])
#    pl.colorbar(ticks=[0.45,0.55])
    pl.savefig('snap'+str(n/skip))
