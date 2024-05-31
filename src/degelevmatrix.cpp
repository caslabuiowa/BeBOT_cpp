#include "../include/degelevmatrix.h"
#include "../include/nchoosek_mod.h"
#include <vector>

std::vector<double> DegElevMatrix(int N, int M) {
    int r = M - N;
    std::vector<std::vector<double>> E(M + 1, std::vector<double>(N + 1, 0));

    for (int i = 1; i <= N + 1; ++i) {
        for (int j = 1; j <= r + 1; ++j) {
            E[i + j - 2][i - 1] = static_cast<double>(nchoosek_mod(N, i - 1)) * 
                                  static_cast<double>(nchoosek_mod(r, j - 1)) / 
                                  static_cast<double>(nchoosek_mod(M, i + j - 2));
        }
    }

    std::vector<double> E_flattened((N + 1) * (M + 1), 0.0);
    for (size_t i = 0; i < E.size(); ++i) {
        for (size_t j = 0; j < E[i].size(); ++j) {
            E_flattened[j * (M + 1) + i] = E[i][j];
        }
    }

    return E_flattened;
}
