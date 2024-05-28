#include "../include/piecewisebernsteinproduct.h"
#include "../include/bernsteinproduct.h"
#include <vector>
#include <algorithm>

std::vector<double> PiecewiseBernsteinProduct(const std::vector<double>& A, const std::vector<double>& B, int K, int N) {
    std::vector<std::vector<double>> A_m(K), B_m(K), Cp_m(K);
    int rows = N + 1;
    for (int i = 0; i < K; ++i) {
        A_m[i] = std::vector<double>(A.begin() + i * rows, A.begin() + (i + 1) * rows);
        B_m[i] = std::vector<double>(B.begin() + i * rows, B.begin() + (i + 1) * rows);
    }

    for (int i = 0; i < K; ++i) {
        Cp_m[i] = BernsteinProduct(A_m[i], B_m[i]);
    }

    std::vector<double> Cp;
    for (auto& cp : Cp_m) {
        Cp.insert(Cp.end(), cp.begin(), cp.end());
    }

    return Cp;
}
