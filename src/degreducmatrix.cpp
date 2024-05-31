#include "../include/degreducmatrix.h"
#include "../include/nchoosek_mod.h"
#include <vector>
#include <Eigen/Dense>  
#include <iostream>

std::vector<std::vector<double>> DegRedMatrix(int M, int N) {
    int r = M - N;
    std::vector<std::vector<double>> E(M + 1, std::vector<double>(N + 1, 0));

    for (int i = 1; i <= N + 1; ++i) {
        for (int j = 1; j <= r + 1; ++j) {
            E[i + j - 2][i - 1] = static_cast<double>(nchoosek_mod(N, i - 1)) * 
                                  static_cast<double>(nchoosek_mod(r, j - 1)) / 
                                  static_cast<double>(nchoosek_mod(M, i + j - 2));
        }
    }

    // Transpose of E
    std::vector<std::vector<double>> E_transposed(N + 1, std::vector<double>(M + 1, 0));
    for (size_t i = 0; i < E.size(); ++i) {
        for (size_t j = 0; j < E[i].size(); ++j) {
            E_transposed[j][i] = E[i][j];
            std::cout << "E_transposed[" << j << "][" << i << "] = " << E_transposed[j][i] << std::endl;

        }
    }

    // Converting to Eigen matrix for pseudo-inverse
    Eigen::MatrixXd Eelev(N + 1, M + 1);
    for (size_t i = 0; i < E_transposed.size(); ++i) {
        for (size_t j = 0; j < E_transposed[i].size(); ++j) {
            Eelev(i, j) = E_transposed[i][j];
            std::cout << "Eelev(" << i << ", " << j << ") = " << Eelev(i, j) << std::endl;
        }
    }

    Eigen::MatrixXd E_pinv = Eelev.completeOrthogonalDecomposition().pseudoInverse();

    // Converting back to std::vector<std::vector<double>>
    std::vector<std::vector<double>> E_result(N + 1, std::vector<double>(M + 1));
    for (int i = 0; i < E_result.size(); ++i) {
        for (int j = 0; j < E_result[i].size(); ++j) {
            E_result[i][j] = E_pinv(i, j);
            std::cout << "E_result[" << i << "][" << j << "] = " << E_result[i][j] << std::endl;
        }
    }

    return E_transposed;
}



// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp$ g++ -I./include -I/usr/include/eigen3 -o myProgram main.cpp bebot/bebot.cpp bebot/bernsteindifferentialmatrix.cpp bebot/bernsteinintegrationmatrix.cpp bebot/bernsteinmatrix_a2b.cpp bebot/bernsteinpoly.cpp bebot/degelevmatrix.cpp bebot/degreducmatrix.cpp bebot/lagrangepoly.cpp bebot/lgl_ps.cpp bebot/nchoosek_mod.cpp bebot/piecewisebebot.cpp bebot/piecewisebernsteinpoly.cpp
