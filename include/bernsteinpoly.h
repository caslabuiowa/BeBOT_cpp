#ifndef BERNSTEINPOLY_H
#define BERNSTEINPOLY_H

#include <vector>

std::vector<std::vector<double>> BernsteinPoly(const std::vector<std::vector<double>>& Cp, const std::vector<double>& t, double t0 = 0.0, double tf = 1.0);

#endif
