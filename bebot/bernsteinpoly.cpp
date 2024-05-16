#include "../include/bernsteinpoly.h"
#include "../include/bernsteinmatrix_a2b.h" 
#include <cmath>
#include <iostream>

template <typename T>
auto BernsteinPoly(const T& Cp, const std::vector<double>& t, double t0, double tf)
    -> typename std::enable_if<std::is_same<T, std::vector<double>>::value, std::vector<double>>::type {
    
    int N = Cp.size() - 1;
    std::vector<double> poly_t(t.size());

    std::vector<std::vector<double>> B = BernsteinMatrix_a2b(N, t, t0, tf);

    for (size_t i = 0; i < t.size(); ++i) {
        poly_t[i] = 0.0;
        for (int k = 0; k < N + 1; ++k) {
            poly_t[i] += Cp[k] * B[i][k];
        }
    }

    return poly_t;
}

template <typename T>
auto BernsteinPoly(const T& Cp, const std::vector<double>& t, double t0, double tf)
    -> typename std::enable_if<std::is_same<T, std::vector<std::vector<double>>>::value, std::vector<std::vector<double>>>::type {
    
    int M = Cp.size();
    int N = Cp[0].size();

    std::vector<std::vector<double>> poly_t(t.size(), std::vector<double>(N));

    // If Cp is a column vector, transpose it
    if (N == 1 && M > 1) {
        N = M;
    }
    N = N - 1;

    std::vector<std::vector<double>> B = BernsteinMatrix_a2b(N, t, t0, tf);

    for (size_t i = 0; i < t.size(); ++i) {
        poly_t[i].resize(1);
        poly_t[i][0] = 0.0;

        for (int k = 0; k < N + 1; ++k) {
            double product = Cp[0][k] * B[i][k];
            poly_t[i][0] += product;
        }
    }

    return poly_t;
}
