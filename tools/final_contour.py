#from numpy import *

import matplotlib
matplotlib.use('Agg')

from pylab import *

def grid(x, y, z , resX=50, resY=50):
    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi , interp='linear')
    X, Y = meshgrid(xi, yi )
    return X, Y, Z

#figure(figsize=( 10.34 , 8))
figure(figsize=( 8 , 8))

#levels = np.linspace( -0.2 , 1.2 , 15 )
#levels = [ -0.2 , -0.15 , - 0.1 , -0.05, 0 , 0.5 , 1 , 1.05 , 1.1 , 1.15, 1.2 ]
levels = [ -0.05, 0 , 0.5 , 1 , 1.05  ]

dtm = loadtxt('200/mesh.dat');
xm=dtm[:,0]; ym=dtm[:,1];  pm=dtm[:,5];  vxm=dtm[:,8]; vym=dtm[:,9]; alm=dtm[:,4]
X, Y, Z = grid(xm, ym, alm)
#contourf( X,Y,Z, levels ,  aspect='equal') #,  colors='g')
CP=contour( X,Y,Z, levels ,
                 colors='k',  # negative contours will be dashed by default
                 aspect='equal') #,  colors='g')

clabel(CP, inline=1, fontsize=10 , fmt='%3.2f') 

xlabel('$x/L$')
ylabel('$y/L$')
xlim([-.7,.7])
ylim([-.7,.7])

#colorbar()

savefig('final_contours',dpi=300)
