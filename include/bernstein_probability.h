#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include "../../BeBOT_cpp/include/bebot.h"
#include "../../BeBOT_cpp/include/bernsteinpoly.h"
#include "../../BeBOT_cpp/include/bernsteinmatrix_a2b.h"

std::vector<std::vector<double>> convertTo2D(const std::vector<double>& flatMatrix,
                                             size_t rows,
                                             size_t cols);

std::vector<std::vector<double>> bernstein_probability(const std::vector<double>& x, 
                                                       const std::vector<double>& x_cdf, 
                                                       const std::vector<double>& P);
