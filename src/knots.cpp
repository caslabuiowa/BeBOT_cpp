#include "../include/knots.h"
#include <vector>

std::vector<double> Knots(const std::vector<double>& Cp, int K) {
    int N = Cp.size() / K - 1;
    std::vector<double> cp_k;
    cp_k.reserve(2 * K);  // Reserve space for 2*K points (first and last of each segment)

    for (int i = 0; i < K; ++i) {
        cp_k.push_back(Cp[i * (N + 1)]);        // First point of each segment
        cp_k.push_back(Cp[(i + 1) * (N + 1) - 1]);  // Last point of each segment
    }

    return cp_k;
}
