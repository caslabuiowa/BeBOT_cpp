#include "../include/piecewiseintegrationmatrix.h"
#include <vector>

std::vector<std::vector<double>> PiecewiseIntegrationMatrix(int K, int N, const std::vector<double>& T) {
    int rows = K * (N + 1);
    int cols = K * (N + 2);
    std::vector<std::vector<double>> I(rows, std::vector<double>(cols, 0));

    for (int i = 0; i < K; ++i) {
        // Calculate the fill value for matrix A
        double fillValue = (T[i + 1] - T[i]) / (N + 1);

        // Fill the values for the triangular part of the matrix I
        for (int j = 1; j < N + 2; ++j) {
            for (int k = 0; k < j; ++k) {
                I[i * (N + 1) + k][(i * (N + 2)) + j] = fillValue;
            }
        }

        // Correctly replicate the MATLAB logic for copying data to the upper right matrix section
        if (i > 0) {
            int startRow = (i - 1) * (N + 1);
            int startCol = i * (N + 2);
            for (int row = startRow; row < startRow + N + 1; ++row) {
                for (int col = startCol; col < cols; ++col) {
                    I[row][col] = I[startRow + N][startCol - 1];
                }
            }
        }
    }
    return I;
}
