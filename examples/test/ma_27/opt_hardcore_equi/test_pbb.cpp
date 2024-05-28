#include "../../../../include/piecewisebebot.h"
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    // Set the order N and the knot vector
    int N = 1;
    std::vector<double> tknots = {0, 2, 4, 6, 8, 10};

    // Create an instance of PiecewiseBeBOT
    PiecewiseBeBOT piecewisebebot(N, tknots);

    // Retrieve the flattened differentiation matrix
    const std::vector<double>& Dm_flat = piecewisebebot.getDifferentiationMatrixFlat();
    std::vector<std::vector<double>> Dm_2D = piecewisebebot.getDifferentiationMatrix();

    // Print the flattened matrix
    std::cout << "Flattened Differentiation Matrix:" << std::endl;
    int size = N + 1;
    int totalSize = (N + 1) * tknots.size() - 1;
    for (int i = 0; i < totalSize; i++) {
        if (i % size == 0 && i != 0 && i % ((N + 1) * 2) == 0) {
            std::cout << std::endl;  // Extra newline to separate segments
        }
        if (i < Dm_flat.size()) {
            std::cout << std::setw(10) << Dm_flat[i] << " ";
        } else {
            std::cout << std::setw(10) << 0 << " ";  // Print zeros for visual consistency
        }
        if ((i + 1) % size == 0) {
            std::cout << std::endl;
        }
    }

    // Print the 2D differentiation matrix
    std::cout << "\n2D Differentiation Matrix:" << std::endl;
    for (int i = 0; i < Dm_2D.size(); i++) {
        for (int j = 0; j < Dm_2D[i].size(); j++) {
            std::cout << std::setw(10) << Dm_2D[i][j] << " ";
        }
        std::cout << std::endl;  // Newline at each row end
    }

    return 0;
}
