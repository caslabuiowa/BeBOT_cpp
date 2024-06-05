#ifndef DE_CASTELJAU_H
#define DE_CASTELJAU_H

#include <vector>
#include <stdexcept>
#include <iostream>

// Type alias for convenience
using Matrix = std::vector<std::vector<double>>;

// Function to create an identity matrix of size n x n
Matrix eye(int n);

// Function to multiply matrix by a scalar
Matrix scalarMultiply(const Matrix& mat, double scalar);

// Function to add two matrices
Matrix addMatrices(const Matrix& A, const Matrix& B);

// Function to multiply two matrices
Matrix multiplyMatrices(const Matrix& A, const Matrix& B);

// Function to extract columns from a matrix
Matrix extractColumns(const Matrix& mat, int colStart, int colEnd);

// deCasteljau function in C++
std::pair<Matrix, std::vector<double>> deCasteljau(const Matrix& Cp, double lambda);

#endif // DE_CASTELJAU_H