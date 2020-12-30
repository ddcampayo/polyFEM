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
    n = str(0)
else:
    n = sys.argv[1]

plt.figure(figsize=(8,8))

path='./'

LL= 1

omega = 2 * np.pi / (31 * 0.01)
print('omega = ' , omega)

Delta_t = 1  #optional prefactor

dt=np.loadtxt(path + n +'/particles.dat')

x=dt[:,0]; y=dt[:,1];
#    vol=dt[:,3]
#w=dt[:,4];
vx=dt[:,8]; vy=dt[:,9];
#p=dt[:,9] / Delta_t**2
#  s=dt[:,10]
#  I=dt[:,11];

#p = 0.5*omega**2 * w

r = np.sqrt( x**2 + y**2 )
v = np.sqrt( vx**2 + vy**2 )

#make furthest pressure value 0

rm = np.argmax(r)

# p -= p[ rm ] #  np.min( p )

plt.plot( r , v , 'o' )

def vel(r) : # analytic solution for Gresho's vortex velocity
    r0 = 0.2

    x = r/r0

    if( x < 1)  : return x
    if( x < 2)  : return 2-x

    return 0;

v_vel = np.vectorize( vel )
rr = np.linspace( 0 , max(r) , 200 )
plt.plot( rr , v_vel(rr)  )

#plt.plot( r , w , 'x' )
   
#plt.xlim([-LL/2.0 , LL/2.0 ])
#plt.ylim([-LL/2.0 , LL/2.0 ])
    #    pl.colorbar(ticks=[0.45,0.55])

#print( 'step no ' + n )

plt.savefig( 'velocity_' + n)
plt.show()

