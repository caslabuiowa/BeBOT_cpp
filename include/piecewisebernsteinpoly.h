#ifndef PIECEWISEBERNSTEINPOLY_H
#define PIECEWISEBERNSTEINPOLY_H

#include <vector>

// Updated function signature to accept a vector of doubles for tknots
std::vector<std::vector<double>> PiecewiseBernsteinPoly(const std::vector<std::vector<double>>& Cp, 
                                                        const std::vector<double>& tknots, 
                                                        const std::vector<double>& t);

#endif
