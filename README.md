# polyFEM

A free and open source pFEM simulation tool, based on the CGAL libs.

Also covered is the procedure of "assignment": transferring information
from a set of moving particles onto a fixed mesh, and back

# Case select For different simulations, simply copy any of the
main_*.cpp files on the _mains_ dir onto main.cpp.

Then, run "source cmake-4.8.sh", and "make"

Sample simu_*.cfg input files are provided, copy them onto simu.cfg in
order to run.

Suggested command to run:

main >& log&


# Features
Currently provided:

Kolmogorov shear flow (main_Kolmo.cpp)
The rotating Zalesak disk test (main_Zalesak.cpp)
Taylor-Green vortex sheet (main_TG.cpp)

# Requires

CGAL 4.8

eigen3 3.2.9 libs:
install as explained in INSTALL file (second option, the one that involves cmake), then
sudo cp -r cmake/ /usr/local/include/eigen3
(On newer versions, it seems simply ```sudo make install``` in the final step works.)

Optional: The suitesparse suite, specifically, cholmod

[![DOI](https://zenodo.org/badge/64474373.svg)](https://zenodo.org/badge/latestdoi/64474373)
