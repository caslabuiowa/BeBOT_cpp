#include "../include/piecewisebebot.h"
#include "../include/bernsteindifferentialmatrix.h"
#include "../include/degelevmatrix.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <iomanip> 

PiecewiseBeBOT::PiecewiseBeBOT(int N, const std::vector<double>& tknots)
    : N(N), originalTknots(tknots) {}


void PiecewiseBeBOT::transformTknots() {
    int M = originalTknots.size() - 1;
    //std::cout << "M: " << M << std::endl;
    transformedTknots.clear();

    for (int i = 0; i < M; ++i) {
        transformedTknots.push_back({originalTknots[i], originalTknots[i + 1]});
    }
}


void PiecewiseBeBOT::calculate() {
    transformTknots();
    int M = transformedTknots.size();
    double T = originalTknots.back() - originalTknots.front();

    tnodes = std::vector<double>((N+1) * M, std::numeric_limits<double>::min());


    // Print tnodes after initialization
    //std::cout << "tnodes after initialization:" << std::endl;
    //for (const auto &node : tnodes) {
    //    std::cout << node << " ";
    //}
    std::cout << std::endl;


    //std::cout << "Contents of transformedTknots:" << std::endl;
    //for (const auto& knotVec : transformedTknots) { // Corrected line
    //    for (const auto& knot : knotVec) {
    //        std::cout << knot << " ";
    //    }
    //    std::cout << std::endl;
    //}

    
    Dm = std::vector<std::vector<double>>((N + 1) * M, std::vector<double>((N + 1) * M, 0.0));

    int blockStartRow = 0;
    int tnodeIndex = 0;
    for (int i = 0; i < M; ++i) {
        double start = transformedTknots[i].front();
        double end = transformedTknots[i].back();

        for (int j = 0; j <= N; ++j, ++tnodeIndex) {
            // Overlapping node included for each segment
            tnodes[tnodeIndex] = start + j * (end - start) / N;
            //std::cout << "tnodes[" << tnodeIndex << "] = " << tnodes[tnodeIndex] << std::endl;
        }

        std::vector<std::vector<double>> Dm_temp = BernsteinDifferentiationMatrix(N, end - start);
        std::vector<std::vector<double>> ElevMatrix = DegElevMatrix(N - 1, N);
        
        // Compute block diagonal matrix
        std::vector<std::vector<double>> Dm_block(N + 1, std::vector<double>(N + 1, 0.0));
        for (size_t i = 0; i < Dm_block.size(); ++i) {
            for (size_t j = 0; j < Dm_block[i].size(); ++j) {
                for (size_t k = 0; k < ElevMatrix.size(); ++k) {
                    if (k < Dm_temp[i].size() and j < ElevMatrix[k].size()) {
                        Dm_block[i][j] += Dm_temp[i][k] * ElevMatrix[k][j];
                    }
                }
            }
        }

        // Place the block in the larger matrix
        for (size_t r = 0; r < Dm_block.size(); ++r) {
            for (size_t c = 0; c < Dm_block[r].size(); ++c) {
                Dm[blockStartRow + r][blockStartRow + c] = Dm_block[r][c];
            }
        }

        // Update the starting row for the next block
        blockStartRow += N + 1;
    }

    w = std::vector<double>((N + 1) * M, T / ((N + 1) * M));
    // Output weights
    //std::cout << "Weights:" << std::endl;
    //for (const auto &weight : w) {
    //    std::cout << std::setprecision(4) << weight << " ";
    //}
    std::cout << "\n\n";

    // Output differentiation matrix
    //std::cout << "Differentiation Matrix:" << std::endl;
    //for (const auto &row : Dm) {
    //    for (const auto &elem : row) {
    //        std::cout << std::setprecision(4) << elem << "\t";
    //    }
    //    std::cout << std::endl;
    //}
    std::cout << "\n";
}

std::vector<double> PiecewiseBeBOT::getNodes() const {
    return tnodes;
}

std::vector<double> PiecewiseBeBOT::getWeights() const {
    return w;
}

std::vector<std::vector<double>> PiecewiseBeBOT::getDifferentiationMatrix() const {
    return Dm;
}
