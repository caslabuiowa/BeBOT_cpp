#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include "bernstein_probability.h"

std::pair<std::vector<double>, std::vector<double>> ecdf(const std::vector<double>& data);
std::tuple<std::pair<double, double>, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> bernstein_estimator(const std::vector<std::pair<double, double>>& estimate);
