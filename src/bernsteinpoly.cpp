#include "../include/bernsteinpoly.h"
#include "../include/bernsteinmatrix_a2b.h" 
#include <cmath>
#include <iostream>

std::vector<std::vector<double>> BernsteinPoly(const std::vector<std::vector<double>>& Cp, const std::vector<double>& t, double t0, double tf) {
    int M = Cp.size();
    int N = Cp[0].size();

    std::vector<std::vector<double>> poly_t(t.size(), std::vector<double>(N));

    // If Cp is a column vector, transpose it
    if (N == 1 && M > 1) {
        N = M;
    }
    N = N - 1;
    //std::cout << "N: " << N << std::endl;
    //std::cout << "M: " << M << std::endl;

    std::vector<std::vector<double>> B = BernsteinMatrix_a2b(N, t, t0, tf);
    // Debug
    //std::cout << "Cp size: " << Cp.size() << "x" << (Cp.empty() ? 0 : Cp[0].size()) << std::endl;
    //std::cout << "B size: " << B.size() << "x" << (B.empty() ? 0 : B[0].size()) << std::endl;
    
    // Print Cp vector
    //std::cout << "Cp vector:" << std::endl;
    //for (const auto& row : Cp) {
    //    for (const auto& elem : row) {
    //        std::cout << elem << " ";
    //    }
    //    std::cout << std::endl;
    //}

    // printing B matrix
    //std::cout << "B matrix:" << std::endl;
    //for (size_t i = 0; i < B.size(); ++i) {
    //    for (size_t j = 0; j < B[i].size(); ++j) {
    //        std::cout << B[i][j] << " ";
    //    }
    //    std::cout << std::endl;
    //}
    
    // Matrix multiplication Cp * B
    for (size_t i = 0; i < t.size(); ++i) { 
    poly_t[i].resize(1);  
    poly_t[i][0] = 0.0;  

    for (int k = 0; k < N+1; ++k) {
        // Accumulate the product
        double product = Cp[0][k] * B[i][k];
        poly_t[i][0] += product;
    }
    // Printing poly_t vector
    //std::cout << "poly_t vector:" << std::endl;
    //for (const auto& row : poly_t) {
    //    for (const auto& elem : row) {
    //        std::cout << elem << " ";
    //    }
    //    std::cout << std::endl;
    //}
}
    
    return poly_t;
}
