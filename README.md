# polyFEM

A free and open source pFEM simulation tool, based on the CGAL libs.

This branch covers the procedure of "assignment": transferring information
from a set of moving particles onto a fixed mesh, and back

# Case select
For different simulations, simply copy any of the main_*.cpp files onto
main.cpp. Same with periodic_*.h onto periodic.h

Then, run "source cmake-4.8.sh", and "make"

Sample simu_*.cfg input files are provided, copy them onto simu.cfg in
order to run.

Suggested command to run:

main >& log&


# Features
Currently provided:

The rotating Zalesak disk test (main_Zalesak.cpp)

# Requires

CGAL 4.8
eigen3 libs


