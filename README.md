# BeBOT_cpp
This repository contains the implementation of the Bernstein/Bezier Optimal Trajectories toolkit [BeBOT](https://github.com/caslabuiowa/BeBOT). 

In the examples you can find implementation of lgl, bebot and piecewisebebot.

## Requirements
To run this examples you must have the following installed:
- CMake (https://cmake.org/)
- Ninja (https://ninja-build.org/)
- Ipopt (https://github.com/coin-or)
- Intel MKL libraries (https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
- HSL libraries (https://licences.stfc.ac.uk/product/coin-hsl-archive)

We reccomed following the following guide to install Ipopt and its dependencies https://coin-or.github.io/Ipopt/INSTALL.html

## Contributors
- [Wladimir Petrov](https://github.com/wladimirpetrov)
- The BeBOT_cpp library is based on the original MATLAB library by [Venanzio Cichella](https://github.com/caslabuiowa/BeBOT_MATLAB) and Python library by [magicbycalvin](https://github.com/caslabuiowa/BeBOT)

Ipopt taken from - [coin-or](https://github.com/coin-or)

MA27, MA57 solvers - (https://www.hsl.rl.ac.uk/catalogue/ma57.html)

Pardiso solver (academic license) - (https://panua.ch/pardiso/)

Spral solver - [ralna](https://github.com/ralna) (SPRAL was built on one CPU only, the next update will include parallel processing).

MKL library - (library https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
