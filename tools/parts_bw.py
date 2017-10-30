#!/usr/bin/python

LL= 2*128.0

import pylab as pl
import sys

pl.figure(figsize=(8,8))

if len(sys.argv) == 1:
   n = 150
else:
   n = int(sys.argv[1])

dtm=pl.loadtxt(str(n)+'/particles.dat')
xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

pl.xlim([ -LL/2 , LL/2 ])
pl.ylim([ -LL/2 , LL/2 ])

size = LL/ pl.sqrt( pl.size(alm) ) * 2.0 #/ 2.0
#line_size = 0.1*size

pl.tick_params(
#    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right='off',
    left='off',
    labelbottom='off',
    labelleft='off')

pl.scatter( xm , ym , c=alm, s = size , linewidths = 0 , cmap='Greys' )

# pl.colorbar()
#    pl.colorbar(ticks=[0.45,0.55])
#    pl.savefig('parts'+str(n/skip))

pl.savefig('parts_bw_'+str(n))

#  pl.scatter( xm , ym , c=alm, s=20*LL/64.0 , linewidths=0 , cmap= pl.get_cmap(name="binary"))
