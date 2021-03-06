#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import sys

#Dt=0.01

#for line in open("simu.cfg"):
# if "step" in line:
#  Dt=float(line.split()[1])

if len(sys.argv) == 1:
   step='0'
else:
   step=sys.argv[1]

#dtx = np.loadtxt(step+'/mesh.dat', dtype = np.float64)
dtx = np.loadtxt('./mesh.dat', dtype = np.float64)

ax=dtx[ : , 27]
ay=dtx[ : , 28]

aa=np.concatenate( ( ax, ay ) )

nbins = int(aa.size/100)

ha_x, bins = np.histogram(aa, nbins, density=1)

h_size = ha_x.size

pdf_aa_x = np.zeros(h_size)
pdf_aa_y = np.zeros(h_size)

outf = open('accel_hist_' + step+ '.dat','w')

for k in range(h_size):
    pdf_aa_x[k] = 0.5*( bins[k]+bins[k+1] )
    pdf_aa_y[k] = ha_x[k]

    outf.write( ' %g %g \n' % (pdf_aa_x[k] , pdf_aa_y[k] ) )

outf.close()


#pl.plot(pdf_ax_x, pdf_ax_y)
#pl.show()

