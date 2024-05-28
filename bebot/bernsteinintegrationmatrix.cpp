#include "../include/bernsteinintegrationmatrix.h"
#include <vector>

std::vector<std::vector<double>> BernsteinIntegrationMatrix(int N, double T) {
    std::vector<std::vector<double>> I(N + 1, std::vector<double>(N + 2, 0));
    double factor = T / (N + 1);

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= i; ++j) {
            I[j][i+1] = factor;
        }
    }

    return I;
}
