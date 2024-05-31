#include "../include/bernsteindifferentialmatrix.h"
#include "mkl.h"
#include <vector>
#include <algorithm>

std::vector<double> BernsteinDifferentiationMatrix(int N, double T) {
    int size = (N + 1) * (N + 1);
    std::vector<double> Dm_flat(size, 0.0);

    for (int i = 0; i < N; ++i) {
        Dm_flat[i * (N + 1) + i] = -N / T; 
        if (i + 1 < N + 1) {
            Dm_flat[(i + 1) * (N + 1) + i] = N / T; 
        }
    }

    return Dm_flat;
}