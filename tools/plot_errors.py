import pylab as pl
from matplotlib import pyplot

dt=pl.loadtxt('timings.dat')

T=dt[:,0];
e=pl.sqrt(dt[:,1]);
dt=dt[:,2];

fig = pyplot.plot()#2,1,1)

pyplot.plot(dt,e, 'o', color='blue')

dt1=pl.array([6e-4,5e-3]);

law1= (2e-3) * (dt1/8e-4)**1;

pyplot.plot(dt1,law1, color='blue', lw=1)

dtp=pl.loadtxt('timings_femP.dat')

Tp=dtp[:,0];
ep=pl.sqrt(dtp[:,1]);
dtp=dtp[:,2];

pyplot.plot(dtp,ep, 'x', color='red')

law2= (4e-2) * (dt1/8e-4)**1;

pyplot.plot(dt1,law2, color='red', lw=1)

dtf=pl.loadtxt('../../quad_again/taylor-green/timing_reg/time_vs_error_t1_fem.dat');

Tf=dtf[:,0];
ef=pl.sqrt(dtf[:,1]);
dtf=dtf[:,2];

pyplot.plot(dtf,ef, '*', color='green')

pyplot.xscale('log')
pyplot.yscale('log')

pyplot.xlabel('$\Delta t$')
pyplot.ylabel('Relative $L_2$ error')

pyplot.show()

pyplot.savefig('fig_err.eps', format='eps', dpi=300)
