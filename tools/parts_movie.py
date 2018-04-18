#!/usr/bin/python

LL=1.0


for line in open("simu.cfg"):
 if "step" in line:
  Dt=float(line.split()[1])
# if "every" in line:
#  skip=int(line.split()[1])

import pylab as pl

pl.figure(figsize=(8,8))

import glob

dirs=glob.glob('[0-9]*')

dirs.sort( key = float )

import pylab as pl

for dir_step in dirs[1:] :

    step = int(dir_step)

    time=Dt*step #+ Dt/2.0
    pl.clf()

    dtm=pl.loadtxt(str(step)+'/mesh.dat')
    #    dtm=pl.loadtxt(path+str(n)+'/mesh.dat')

    xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];
    vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

    vel=vxm**2+vym**2

    #    pl.scatter( xm , ym , c=alm, s=80, vmin= -limits , vmax= limits)
    pl.scatter( xm , ym , c= pm , s=80/LL , cmap= pl.cm.bwr , linewidths=0 )
    pl.colorbar()
    pl.quiver( xm , ym , vxm , vym )

    pl.xlim([ -LL/2 , LL/2 ])
    pl.ylim([ -LL/2 , LL/2 ])
    #    pl.colorbar(ticks=[0.45,0.55])
    #    pl.savefig('parts'+str(n/skip)  )

    pl.savefig('parts'+str(step).zfill(5) )

#  pl.scatter( xm , ym , c=alm, s=20*LL/64.0 , linewidths=0 , cmap= pl.get_cmap(name="binary"))

