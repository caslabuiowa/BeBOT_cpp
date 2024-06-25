#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>

// Include the necessary headers for your Bernstein probability functions
#include "../include/bernstein_probability.h"
#include "../include/bernstein_estimator.h"

// Function to calculate ECDF (Empirical Cumulative Distribution Function)

std::pair<std::vector<double>, std::vector<double>> ecdf(const std::vector<double>& data) {
    // Sort the data
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());

    // Compute the ECDF values
    std::vector<double> cdf(sorted_data.size() + 1);
    cdf[0] = 0.0; // First element is 0
    for (size_t i = 1; i < sorted_data.size(); ++i) {
        cdf[i] = static_cast<double>(i) / sorted_data.size();
    }
    cdf.back() = 1.0; // Last element is 1

    // Include the first element in sorted_data to match the cdf size
    sorted_data.insert(sorted_data.begin(), sorted_data.front());

    return {cdf, sorted_data};
}

// Main function: BernsteinEstimator
std::tuple<std::pair<double, double>, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>>
bernstein_estimator(const std::vector<std::pair<double, double>>& estimate) {

    std::vector<double> x_meas, y_meas;

    // Separate into x and y vectors
    for (const auto& pair : estimate) {
        x_meas.push_back(pair.first);
        y_meas.push_back(pair.second);
    }

    // Define interval bounds
    const double a_x = -50.0, a_y = -50.0, b_x = 50.0, b_y = 50.0;

    // Project data
    std::vector<double> x_est, y_est;
    std::transform(x_meas.begin(), x_meas.end(), std::back_inserter(x_est), [a_x, b_x](double val) {
        return (val - a_x) / (b_x - a_x);
    });
    std::transform(y_meas.begin(), y_meas.end(), std::back_inserter(y_est), [a_y, b_y](double val) {
        return (val - a_y) / (b_y - a_y);
    });

    // Build estimator
    auto [Fx, x_cdf] = ecdf(x_est);
    auto [Fy, y_cdf] = ecdf(y_est);

    // Compute optimal order for BP
    size_t n = x_meas.size();
    size_t m = static_cast<size_t>(std::ceil(std::pow(n, 4.0 / 5.0))) + 2;

    // Create BP to approximate CDF
    std::vector<double> x(m), y(m);
    std::generate(x.begin(), x.end(), [m, idx = 0]() mutable {
        return static_cast<double>(idx++) / (m - 1);
    });
    std::generate(y.begin(), y.end(), [m, idx = 0]() mutable {
        return static_cast<double>(idx++) / (m - 1);
    });

    auto temp_x = bernstein_probability(x, x_cdf, Fx);
    auto temp_y = bernstein_probability(y, y_cdf, Fy);

    // Rescale
    std::vector<double> rescaled_x(m), rescaled_y(m);
    std::transform(x.begin(), x.end(), rescaled_x.begin(), [a_x, b_x](double val) {
        return val * (b_x - a_x) + a_x;
    });
    std::transform(y.begin(), y.end(), rescaled_y.begin(), [a_y, b_y](double val) {
        return val * (b_y - a_y) + a_y;
    });

    std::vector<double> pdf_x(m), pdf_y(m);
    for (size_t i = 0; i < m; ++i) {
        pdf_x[i] = temp_x[i][0];
        pdf_y[i] = temp_y[i][0];
    }

    // Normalize pdf_x and pdf_y
    std::vector<double> rescaled_pdf_x(m), rescaled_pdf_y(m);
    double divisor_x = b_x - a_x;
    double divisor_y = b_y - a_y;
    std::transform(pdf_x.begin(), pdf_x.end(), rescaled_pdf_x.begin(), [divisor_x](double x) {
        return x / divisor_x;
    });
    std::transform(pdf_y.begin(), pdf_y.end(), rescaled_pdf_y.begin(), [divisor_y](double y) {
        return y / divisor_y;
    });

    // Find average of distribution
    auto max_x_it = std::distance(rescaled_pdf_x.begin(), std::max_element(rescaled_pdf_x.begin(), rescaled_pdf_x.end()));
    auto max_y_it = std::distance(rescaled_pdf_y.begin(), std::max_element(rescaled_pdf_y.begin(), rescaled_pdf_y.end()));
    double new_estimate_x = rescaled_x[max_x_it];
    double new_estimate_y = rescaled_y[max_y_it];

    std::pair<double, double> new_estimate = {new_estimate_x, new_estimate_y};
    std::vector<std::pair<double, double>> prob_x(m), prob_y(m);
    for (size_t i = 0; i < m; ++i) {
        prob_x[i] = std::make_pair(rescaled_x[i], rescaled_pdf_x[i]);
        prob_y[i] = std::make_pair(rescaled_y[i], rescaled_pdf_y[i]);
    }

    return {new_estimate, prob_x, prob_y};
}