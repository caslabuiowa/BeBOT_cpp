#include <iostream>
#include <vector>
#include <stdexcept>
#include <tuple>
#include "/opt/intel/oneapi/mkl/latest/include/mkl.h"

// Type alias for convenience
using Matrix = std::vector<std::vector<double>>;

// Function to create an identity matrix of size n x n
Matrix eye(int n) {
    Matrix I(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        I[i][i] = 1.0;
    }
    return I;
}

// Function to multiply matrix by a scalar
Matrix scalarMultiply(const Matrix& mat, double scalar) {
    Matrix result = mat;
    for (auto& row : result) {
        for (auto& elem : row) {
            elem *= scalar;
        }
    }
    return result;
}

// Function to add two matrices
Matrix addMatrices(const Matrix& A, const Matrix& B) {
    int rows = A.size();
    int cols = A[0].size();
    Matrix C(rows, std::vector<double>(cols));

    for (int i = 0; i < rows; ++i) {
        // Add rows using cblas_daxpy
        cblas_daxpy(cols, 1.0, A[i].data(), 1, C[i].data(), 1);
        cblas_daxpy(cols, 1.0, B[i].data(), 1, C[i].data(), 1);
    }

    return C;
}

// Function to multiply two matrices
Matrix multiplyMatrices(const Matrix& A, const Matrix& B) {
    if (A[0].size() != B.size()) {
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication.");
    }
    Matrix result(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < B[0].size(); ++j) {
            for (size_t k = 0; k < B.size(); ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Function to extract columns from a matrix
Matrix extractColumns(const Matrix& mat, int colStart, int colEnd) {
    Matrix result(mat.size(), std::vector<double>(colEnd - colStart));
    for (size_t i = 0; i < mat.size(); ++i) {
        for (int j = colStart; j < colEnd; ++j) {
            result[i][j - colStart] = mat[i][j];
        }
    }
    return result;
}

Matrix selectSubmatrix(const Matrix& mat, const std::vector<int>& row_indices, const std::vector<int>& col_indices) {

    int row_start = row_indices.front();
    int row_end = row_indices.back();

    int col_start = col_indices.front();
    int col_end = col_indices.back();

    if (row_start < 0 || row_end > mat.size() || col_start < 0 || col_end > mat[0].size()) {
        throw std::out_of_range("Index out of bounds");
    }

    int num_rows = row_end - row_start + 1;
    int num_cols = col_end - col_start + 1;
    Matrix submatrix(num_rows, std::vector<double>(num_cols));

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            submatrix[i][j] = mat[row_start + i][col_start + j];
        }
    }

    return submatrix;
}

void replaceSubmatrix(Matrix& mat, const Matrix& submat, size_t row_start, size_t col_start) {
    size_t num_rows = submat.size();
    size_t num_cols = submat[0].size();

    // Ensure the submatrix fits within the larger matrix
    if (row_start + num_rows > mat.size() || col_start + num_cols > mat[0].size()) {
        throw std::out_of_range("Submatrix exceeds bounds of the larger matrix");
    }

    // Replace the submatrix
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            mat[row_start + i][col_start + j] = submat[i][j];
        }
    }
}

// deCasteljau function in C++
std::tuple<Matrix, Matrix, std::vector<double>> deCasteljau(const Matrix& Cp, double lambda) {
    size_t m = Cp.size();
    size_t n = Cp[0].size() - 1;

    Matrix A = scalarMultiply(eye(n), 1 - lambda);
    Matrix B = scalarMultiply(eye(n), lambda);
    Matrix C(n + 1, std::vector<double>(n, 0.0));
    

    for (size_t row = 0; row < n; ++row) {
        for (size_t col = 0; col < n; ++col) {
            C[row][col] = C[row][col] + A[row][col];
        }
    }

    for (size_t row = 0; row < n; ++row) {
        for (size_t col = 0; col < n; ++col) {
            C[row + 1][col] = C[row + 1][col] + A[row][col];
        }
    }

    for (const auto& row : C) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    Matrix Cpout(m, std::vector<double>(2 * n + 1));
    Matrix Cptemp = Cp;
    Matrix Cp1(m, std::vector<double>(n));
    Matrix Cp2(m, std::vector<double>(n));

    std::vector<int> row_indices = {0, 1};
    std::vector<int> col_indices = {1, 3};

    // Select the specified rows and columns
    Matrix Cptemp_new = selectSubmatrix(Cptemp, row_indices, col_indices);
    for (const auto& row : Cptemp) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    for (const auto& row : Cptemp_new) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    for (int i = n; i >= 1; --i) {
        Matrix temp = multiplyMatrices(selectSubmatrix(Cptemp,{0, Cptemp.size()-1},{0, i}), selectSubmatrix(C,{0, i},{0, i-1}));
        replaceSubmatrix(Cptemp, temp, 0, 0);
        for (int j = 0; j < m; ++j) {
            Cp1[j][n - i] = Cptemp[j][0];
            Cp2[j][i - 1] = Cptemp[j][i - 1];
        }
    }

    // Assign first and last control points
    for (int i = 0; i < m; ++i) {
        Cp1[i].insert(Cp1[i].begin(),Cp[i][0]);
        Cp2[i].push_back(Cp[i][n]);
    }

    std::vector<double> Pos(m);
    for (int i = 0; i < m; ++i) {
        Pos[i] = Cp1[i][n - 1];
    }

    return std::make_tuple(Cp1, Cp2, Pos);
}