#include "../include/bebot.h"
#include "../include/bernsteindifferentialmatrix.h"
#include "../include/degelevmatrix.h"
#include <iostream>
#include <iomanip>
#include "mkl.h"


Bebot::Bebot(int N, double T) : N(N), T(T) {}

void Bebot::calculate() {
    // tnodes calculation
    tnodes.resize(N + 1);
    double step = T / N;
    for (int i = 0; i <= N; ++i) {
        tnodes[i] = i * step;
    }

    // w calculation
    w.resize(N + 1, T / (N + 1));

    // dm calculation
    //auto Dm_temp = BernsteinDifferentiationMatrix(N, T);
    std::vector<double> Dm_temp_flat = BernsteinDifferentiationMatrix(N, T);
    //auto ElevMatrix = DegElevMatrix(N - 1, N);
    std::vector<double> ElevMatrix_flat = DegElevMatrix(N - 1, N);
    
    // Print Dm_temp
    //std::cout << "\nDm_temp Matrix:" << std::endl;
    //for (const auto &row : Dm_temp) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setprecision(4) << elem << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    // Flatten Dm_temp
    //std::vector<double> Dm_temp_flat((N + 1) * (N + 1));
    //for (size_t i = 0; i < Dm_temp.size(); ++i) {
    //    for (size_t j = 0; j < Dm_temp[i].size(); ++j) {
    //        Dm_temp_flat[j + i * (N + 1)] = Dm_temp[i][j];
    //    }
    //}

    // Flatten ElevMatrix
    //std::vector<double> ElevMatrix_flat((N + 1) * (N + 1));
    //for (size_t i = 0; i < ElevMatrix.size(); ++i) {
    //    for (size_t j = 0; j < ElevMatrix[i].size(); ++j) {
    //        ElevMatrix_flat[j + i * (N + 1)] = ElevMatrix[i][j];
    //    }
    //}
    
    
    // Initializing size of Dm 
    Dm.resize((N + 1) * (N + 1), 0.0);
    

    // Perform matrix multiplication using MKL
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 
                N + 1, N + 1, N + 1, 
                1.0, 
                &Dm_temp_flat[0], N + 1, 
                &ElevMatrix_flat[0], N + 1, 
                0.0, 
                &Dm[0], N + 1);
    // Printing Dm
    //std::cout << "\nDm Matrix:" << std::endl;
    //for (const auto &row : Dm) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setprecision(4) << elem << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    // Initializing Dm with the correct size and zero values
    //Dm = std::vector<std::vector<double>>(N + 1, std::vector<double>(N + 1, 0));

    // Dm = Dm_temp * ElevMatrix
    //for (size_t i = 0; i <= N; ++i) {
    //    for (size_t j = 0; j <= N; ++j) {
    //        for (size_t k = 0; k <= N; ++k) {
    //            if (k < Dm_temp[i].size() && j < ElevMatrix[k].size()) {
    //                Dm[i * (N + 1) + j] += Dm_temp[i][k] * ElevMatrix[k][j];
    //            }
    //        }
    //    }
    //}

    // Printing Dm after it has been fully calculated
    //std::cout << "\nFinal Dm Matrix:" << std::endl;
    //for (const auto &row : Dm) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setw(10) << std::setprecision(4) << elem << " ";
    //    }
    //    std::cout << std::endl;
    //}
}

std::vector<double> Bebot::getNodes() {
    return tnodes;
}

std::vector<double> Bebot::getWeights() {
    return w;
}

const std::vector<double>& Bebot::getDifferentiationMatrix() const {
    return Dm;
}
