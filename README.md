# BeBOT_cpp
This repository contains the implementation of the Bernstein/Bezier Optimal Trajectories toolkit [BeBOT](https://github.com/caslabuiowa/BeBOT). Python and Matlab implementations of this library are also available here [Python BeBOT](https://github.com/caslabuiowa/BeBOT) and [MATLAB BeBOT](https://github.com/caslabuiowa/BeBOT_MATLAB).

In the examples you can find implementation of lgl, bebot and piecewisebebot.

## Requirements
To run this examples you must have the following installed:
- CMake (https://cmake.org/)
- Ninja (https://ninja-build.org/)
- Ipopt (https://github.com/coin-or)
- Intel MKL libraries (https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
- HSL libraries (https://licences.stfc.ac.uk/product/coin-hsl-archive)

We reccomed following the following guide to install Ipopt and its dependencies https://coin-or.github.io/Ipopt/INSTALL.html

## Example 1
In the example provided, the following optimal control problem is solved using Bernstein polynomials.
$$
\min_{x(t),u(t),t_f} \int_0^{t_f} 1 d\tau
$$