#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize

import sys

#print "This is the name of the script: ", sys.argv[0]
#print "Number of arguments: ", len(sys.argv)
#print "The arguments are: " , str(sys.argv)

if(len(sys.argv) == 1) :
    init_t = 0
else:
    init_t = int( sys.argv[1] )

#import pylab as pl


plt.figure(figsize=(8,8))

skip=1
#path='timings_full/'
path='./'

LL= 1

for n in range( init_t ,2000000+skip,skip):
    plt.clf()
    dt=np.loadtxt(path+str(n)+'/particles.dat')
 
    x=dt[:,0]; y=dt[:,1];
    vol=dt[:,3]
#    vx=dt[:,5]; vym=dt[:,6];
    p=dt[:,5]

#    I=dt[:,14];  #  eccentricity



    r = np.sqrt( x**2 + y**2 )

    rm = np.argmax(r)

    p -= p[ rm ] #  np.min( p )



    #    plt.plot( r , p , 'o' )

    plt.scatter( x , y , s=80 , c=p )
#    plt.scatter( x , y , 80, c= vol , vmin=0.0022, vmax=0.0028 )
#    plt.scatter( x , y , 10, c=w )
#    plt.scatter( x , y , 10, c=I )
#    plt.scatter( x , y , 80,  c= I , vmin= 1.02e-6, vmax= 1.06e-6 )
#    plt.scatter( x , y , 80,  c= np.log( d2 + 1e-18 ) )
#    plt.scatter( x , y , 10, c=om )

    plt.xlim([-LL/2.0 , LL/2.0 ])
    plt.ylim([-LL/2.0 , LL/2.0 ])
    #    pl.colorbar(ticks=[0.45,0.55])

    print( 'snap{:03d}'.format( int(n/skip)  ) )

    plt.savefig( 'snap{:03d}'.format( int(n/skip)  ) )

