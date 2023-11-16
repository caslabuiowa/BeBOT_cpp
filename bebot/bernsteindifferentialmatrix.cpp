#include "../include/bernsteindifferentialmatrix.h"
#include <vector>

std::vector<std::vector<double>> BernsteinDifferentiationMatrix(int N, double T) {
    std::vector<std::vector<double>> Dm(N + 1, std::vector<double>(N, 0));

    for (int i = 0; i < N; ++i) {
        Dm[i][i] = -N / T;
        Dm[i + 1][i] = N / T;
    }

    return Dm;
}
