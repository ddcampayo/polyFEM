#!/usr/bin/python

LL= 2*64.0

import pylab as pl
import sys

pl.figure(figsize=(8,8))

pl.tight_layout()

if len(sys.argv) == 1:
   n = 150
else:
   n = int(sys.argv[1])

ax=pl.axes()

dtm=pl.loadtxt(str(n)+'/particles.dat')
xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]

cp = 2   # copies

pl.xlim([ -LL/2 , LL*( cp - 1 + 1.0/2.0) ])
pl.ylim([ -LL/2 , LL*( cp - 1 + 1.0/2.0) ])

#size = 8*LL/64.0
size=10

ax.xaxis.set_major_locator(pl.NullLocator())
ax.yaxis.set_major_locator(pl.NullLocator())

for dx in range(0,cp) :
   for dy in range(0,cp) :

      ax.scatter( xm + dx * LL , ym  + dy * LL , c=alm, s = size , linewidths=0 , cmap='Greys' )

#pl.colorbar()
#    pl.colorbar(ticks=[0.45,0.55])
#    pl.savefig('parts'+str(n/skip))

#pl.savefig( 'parts_bw_values_' + str(n) + '.pdf' , format='pdf' , bbox_inches = 'tight' )

pl.savefig( 'parts_bw_values_4_' + str(n) + '.png' , format='png' , bbox_inches = 'tight' )

#  pl.scatter( xm , ym , c=alm, s=20*LL/64.0 , linewidths=0 , cmap= pl.get_cmap(name="binary"))
