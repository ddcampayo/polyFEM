#!/usr/bin/python
import glob

dirs=glob.glob('./quad?')
#dirs=glob.glob('./femP?')
#dirs=glob.glob('./fem?')
#dirs=glob.glob('./quad_m6_?')

#dirs = sorted(dirs)

dirs.sort( ) # key = char )

for case in dirs :

# print case

 for line in open(case+"/simu.cfg"):
  if "step" in line:
   Dt=float(line.split()[1])

 with open(case+"/log", 'rb') as lf:
  last = lf.readlines()[-1].decode()

 #with open(case+"/L2.dat", 'rb') as l2:
 with open(case+"/L2_p.dat", 'rb') as l2:
  l22 = l2.readlines()[-1].decode()


 print Dt, l22.split()[1] ,   last.split()[2]
# print " %g  %g " % (
#  Dt ,
#  last
# )

