#include "../include/bernstein_probability.h"
#include <vector>
#include <algorithm>
#include <cblas.h>

// Helper function to convert a 1D vector to a 2D vector
std::vector<std::vector<double>> convertTo2D(const std::vector<double>& flatMatrix, size_t rows, size_t cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix[i][j] = flatMatrix[i * cols + j];
        }
    }
    return matrix;
}

std::vector<std::vector<double>> bernstein_probability(const std::vector<double>& x, 
                                                       const std::vector<double>& x_cdf, 
                                                       const std::vector<double>& P) {
    size_t n = x.size();
    std::vector<double> prob(n, 0.0);

    for (size_t k = 0; k < n; ++k) {
        // Find the last index in x_cdf where x_cdf[index] <= x[k]
        auto it = std::find_if(x_cdf.rbegin(), x_cdf.rend(), [&](double val) { return val <= x[k]; });

        if (it == x_cdf.rend()) {
            prob[k] = P[0];
        } else {
            size_t temp_x = std::distance(it, x_cdf.rend()) - 1;
            prob[k] = P[temp_x];
        }
    }
    
    // Assuming BeBOT class is defined and has methods `calculate()` and `getDifferentiationMatrix()`
    Bebot bebot(n-1, 1);
    bebot.calculate();
    const auto& Dm = bebot.getDifferentiationMatrix();

    // Output vector
    std::vector<double> pdf(n, 0.0);
    cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, Dm.data(), n, prob.data(), 1, 0.0, pdf.data(), 1);

    // Assuming BernsteinPoly function is defined
    std::vector<std::vector<double>> pdf_2d(1, pdf);
    std::vector<std::vector<double>> BP_dot = BernsteinPoly(pdf_2d, x, 0, 1);
    return BP_dot;
}