#include "../include/bebot.h"
#include "../include/bernsteindifferentialmatrix.h"
#include "../include/degelevmatrix.h"
#include <iostream>
#include <iomanip>

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
    auto Dm_temp = BernsteinDifferentiationMatrix(N, T);
    auto ElevMatrix = DegElevMatrix(N - 1, N);
    
    
    // Print Dm_temp
    //std::cout << "\nDm_temp Matrix:" << std::endl;
    //for (const auto &row : Dm_temp) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setprecision(4) << elem << "\t";
    //    }
    //    std::cout << std::endl;
    //}
    
    
    // Initializing size of Dm 
    Dm.resize(N + 1, std::vector<double>(N + 1, 0.0));
    
    // Printing Dm
    //std::cout << "\nDm Matrix:" << std::endl;
    //for (const auto &row : Dm) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setprecision(4) << elem << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    // Initializing Dm with the correct size and zero values
    Dm = std::vector<std::vector<double>>(N + 1, std::vector<double>(N + 1, 0));

    // Dm = Dm_temp * ElevMatrix
    for (size_t i = 0; i < Dm.size(); ++i) {
        for (size_t j = 0; j < Dm[i].size(); ++j) {
            for (size_t k = 0; k < ElevMatrix.size(); ++k) {
                if (k < Dm_temp[i].size() && j < ElevMatrix[k].size()) {
                    Dm[i][j] += Dm_temp[i][k] * ElevMatrix[k][j];
                }
            }
        }
    }
}

std::vector<double> Bebot::getNodes() {
    return tnodes;
}

std::vector<double> Bebot::getWeights() {
    return w;
}

std::vector<std::vector<double>> Bebot::getDifferentiationMatrix() {
    return Dm;
}
