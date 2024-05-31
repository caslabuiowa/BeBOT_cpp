#include "../include/bernsteinmatrix_a2b.h"
#include <cmath>
#include <iostream>
#include <iomanip>  

std::vector<std::vector<double>> BernsteinMatrix_a2b(int N, const std::vector<double>& t, double tmin, double tmax) {
    //std::cout << "tmax: " << tmax << std::endl;
    //std::cout << "tmin: " << tmin << std::endl;
    if (!t.empty() && tmin == 0) {
        tmax = tmax;
        //std::cout << "tmax: " << tmax << std::endl;
        tmin = tmin;
        //std::cout << "tmin: " << tmin << std::endl;
    } else {
        // tmax and tmin retain their original values
        //std::cout << "Using provided tmax: " << tmax << " and tmin: " << tmin << std::endl;
    }

    std::vector<double> B(N + 1, 1.0); 
    for (int j = 1; j <= N / 2; ++j) {
        B[j] = B[j - 1] * (N + 1 - j) / j;
        B[N - j] = B[j];
    }

    std::vector<std::vector<double>> b(t.size(), std::vector<double>(N + 1, 1.0));
    for (size_t i = 0; i < t.size(); ++i) {
        double tp = t[i] - tmin;
        double ttp = tmax - t[i];
        std::vector<double> T(N + 1, 1.0);
        std::vector<double> TT(N + 1, 1.0);
        for (int j = 1; j <= N; ++j) {
            T[j] = tp * T[j - 1];
            TT[N - j] = ttp * TT[N - j + 1];
            T[j - 1] *= B[j - 1];
        }
        for (int j = 0; j <= N; ++j) {
            b[i][j] = T[j] * TT[j] / pow(tmax - tmin, N);
        }
    }

    // Printing the matrix b
    //std::cout << "\nBernstein Matrix (b):" << std::endl;
    //for (const auto &row : b) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setprecision(4) << elem << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    return b;
}
