
#include <vector>
#include <iostream>

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
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        throw std::invalid_argument("Matrices must have the same dimensions for addition.");
    }
    Matrix result = A;
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] += B[i][j];
        }
    }
    return result;
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

// deCasteljau function in C++
std::pair<Matrix, std::vector<double>> deCasteljau(const Matrix& Cp, double lambda) {
    int m = Cp.size();
    int n = Cp[0].size() - 1;

    Matrix A = scalarMultiply(eye(n), 1 - lambda);
    Matrix B = scalarMultiply(eye(n), lambda);
    Matrix C(n + 1, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = A[i][j];
            if (i + 1 < n + 1) {
                C[i + 1][j] += B[i][j];
            }
        }
    }

    Matrix Cpout(m, std::vector<double>(2 * n + 1));
    Matrix Cptemp = Cp;
    Matrix Cp1(m, std::vector<double>(n));
    Matrix Cp2(m, std::vector<double>(n));

    for (int i = n; i >= 1; --i) {
        Matrix temp = multiplyMatrices(extractColumns(Cptemp, 0, i + 1), extractColumns(C, 0, i));
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < i; ++k) {
                Cptemp[j][k] = temp[j][k];
            }
        }
        for (int j = 0; j < m; ++j) {
            Cp1[j][n - i] = Cptemp[j][0];
            Cp2[j][i - 1] = Cptemp[j][i - 1];
        }
    }

    for (int i = 0; i < m; ++i) {
        Cpout[i][0] = Cp[i][0];
        for (int j = 1; j < n + 1; ++j) {
            Cpout[i][j] = Cp1[i][j - 1];
        }
        for (int j = n + 1; j < 2 * n; ++j) {
            Cpout[i][j] = Cp2[i][j - n];
        }
        Cpout[i][2 * n] = Cp[i][n];
    }

    std::vector<double> Pos(m);
    for (int i = 0; i < m; ++i) {
        Pos[i] = Cp1[i][n - 1];
    }

    return {Cpout, Pos};
}

// Helper function to print a matrix
void printMatrix(const Matrix& mat) {
    for (const auto& row : mat) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

// Helper function to print a vector
void printVector(const std::vector<double>& vec) {
    for (const auto& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

int main() {
    // Example control points
    Matrix Cp = {
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 0.0, -1.0}
    };

    double lambda = 0.5; // Example lambda value

    // Call deCasteljau function
    auto [Cpout, Pos] = deCasteljau(Cp, lambda);

    // Print the results
    std::cout << "Cpout:" << std::endl;
    printMatrix(Cpout);

    std::cout << "Pos:" << std::endl;
    printVector(Pos);

    return 0;
}