#include "../include/bernsteinproduct.h"
#include "../include/nchoosek_mod.h"
#include <vector>

std::vector<double> BernsteinProduct(const std::vector<double>& A, const std::vector<double>& B) {
    int M = A.size() - 1;
    int N = B.size() - 1;
    std::vector<double> C(M + N + 1, 0);

    for (int k = 0; k <= M + N; ++k) {
        for (int j = std::max(0, k - N); j <= std::min(M, k); ++j) {
            C[k] += nchoosek_mod(M, j) * nchoosek_mod(N, k - j) / static_cast<double>(nchoosek_mod(M + N, k)) * A[j] * B[k - j];
        }
    }

    return C;
}
