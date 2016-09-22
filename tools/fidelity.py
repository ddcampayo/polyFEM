#!/usr/bin/python

import pylab as pl

skip=20
#path='timings_full/'
path='./'



for n in range(0,200+skip,skip):

    dt0=pl.loadtxt('../master/'+str(n)+'/particles.dat')

    al0=dt0[:,4]

    dtm=pl.loadtxt(path+str(n)+'/particles.dat')

#    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9];
    alm=dtm[:,4]
    
    print n , pl.norm( al0 - alm ) / pl.norm(al0 )
