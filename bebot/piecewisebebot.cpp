#include "../include/piecewisebebot.h"
#include "../include/bernsteindifferentialmatrix.h"
#include "../include/degelevmatrix.h"
#include "mkl.h"
#include <iostream>
#include <vector>
#include <cmath>

PiecewiseBeBOT::PiecewiseBeBOT(int N, const std::vector<double>& tknots)
    : N(N), originalTknots(tknots) {
    calculate();
}

void PiecewiseBeBOT::transformTknots() {
    transformedTknots.clear();
    for (size_t i = 0; i < originalTknots.size() - 1; ++i) {
        transformedTknots.push_back({originalTknots[i], originalTknots[i + 1]});
    }
}

void PiecewiseBeBOT::calculate() {
    transformTknots();
    double T = originalTknots.back() - originalTknots.front();
    int M = transformedTknots.size();
    tnodes.resize((N + 1) * M);
    w.resize((N + 1) * M, T / ((N + 1) * M));
    Dm_flat.resize((N + 1) * (N + 1) * M);

    for (int segment = 0; segment < M; ++segment) {
        double start = transformedTknots[segment][0];
        double end = transformedTknots[segment][1];
        double segmentLength = end - start;

        for (int i = 0; i <= N; ++i) {
            tnodes[segment * (N + 1) + i] = start + (segmentLength * i / N);
        }

        std::vector<double> Dm_temp = BernsteinDifferentiationMatrix(N, segmentLength);
        std::vector<double> ElevMatrix = DegElevMatrix(N - 1, N);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N + 1, N + 1, N + 1,
                    1.0, Dm_temp.data(), N + 1,
                    ElevMatrix.data(), N + 1,
                    0.0, &Dm_flat[segment * (N + 1) * (N + 1)], N + 1);
    }
}

std::vector<double> PiecewiseBeBOT::getNodes() const {
    return tnodes;
}

std::vector<double> PiecewiseBeBOT::getWeights() const {
    return w;
}

std::vector<std::vector<double>> PiecewiseBeBOT::getDifferentiationMatrix() const {
    // Convert Dm_flat back to 2D matrix if needed for compatibility
    int totalSize = (N + 1) * originalTknots.size() - 1; // Adjust for segments
    std::vector<std::vector<double>> Dm_2D(totalSize, std::vector<double>(totalSize, 0.0));
    int M = transformedTknots.size();

    for (int segment = 0; segment < M; ++segment) {
        for (int i = 0; i < N + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                Dm_2D[segment * (N + 1) + i][segment * (N + 1) + j] = Dm_flat[segment * (N + 1) * (N + 1) + i * (N + 1) + j];
            }
        }
    }

    return Dm_2D;
}

const std::vector<double>& PiecewiseBeBOT::getDifferentiationMatrixFlat() const {
    return Dm_flat;
}
