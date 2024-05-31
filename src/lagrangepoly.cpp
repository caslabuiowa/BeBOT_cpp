#include "../include/lagrangepoly.h"

std::vector<double> LagrangePoly(const std::vector<double>& x, const std::vector<double>& tnodes, const std::vector<double>& time) {
    size_t n = x.size();
    std::vector<double> xN(time.size(), 0.0);

    for (size_t i = 0; i < n; ++i) {
        std::vector<double> basis_polynomial(time.size(), 1.0);

        for (size_t k = 0; k < n; ++k) {
            if (i != k) {
                for (size_t t = 0; t < time.size(); ++t) {
                    basis_polynomial[t] *= (time[t] - tnodes[k]) / (tnodes[i] - tnodes[k]);
                }
            }
        }

        for (size_t t = 0; t < time.size(); ++t) {
            xN[t] += x[i] * basis_polynomial[t];
        }
    }

    return xN;
}
