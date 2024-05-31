#include "../include/lgl_ps.h"
#include <cmath>
#include <iostream>

LGL_PS::LGL_PS(int N, double tf) : N(N), tf(tf) {}

void LGL_PS::calculate() {
    // Creating vectors for storing the LGL nodes and Legendre polynomials
    tnodes.resize(N + 1);
    w.resize(N + 1);
    Dm.resize(N + 1, std::vector<double>(N + 1));

    // Initialize z values
    std::vector<double> z(N + 1, 0.0);
    for (int i = 0; i <= N; i++) {
        z[i] = cos(M_PI * (i * 2.0 + 1) / (2.0 * N + 2.0));
    }

    std::vector<std::vector<double>> L(N + 1, std::vector<double>(N + 1, 0.0));

    /////////////////////////////////////////////////////
    // Finding roots of (1-z^2)LNdot using Newton-Raphson
    /////////////////////////////////////////////////////

    // Setting the convergence criteria for the Newton-Raphson method
    double eps = 2.2204e-16; // choosing convergence criterion to be a very small value
    double max_diff = 2 * eps;

    // Starting the Newton-Raphson loop to find LGL nodes
    while (max_diff > eps) {
        std::vector<double> zprev = z;

        // Initializing the first two Legendre polynomials
        for (int i = 0; i <= N; i++) {
            L[i][0] = 1.0;
            L[i][1] = z[i];
        }

        // Computing Legendre polynomials using a recurrence relation
        for (int i = 2; i <= N; i++) {
            for (int j = 0; j <= N; j++) {
                L[j][i] = ((2 * i - 1) * z[j] * L[j][i - 1] - (i - 1) * L[j][i - 2]) / i;
            }
        }

        // Computing the Newton-Raphson update for z
        for (int i = 0; i <= N; i++) {
            double num = (N + 1) * (z[i] * L[i][N] - L[i][N - 1]);
            double den = (N + 1) * (N + 2) * L[i][N];
            z[i] = zprev[i] - num / den;
        }

        // Calculating the maximum difference between current and previous z values
        max_diff = 0.0;
        for (int i = 0; i <= N; i++) {
            max_diff = std::max(max_diff, std::abs(z[i] - zprev[i]));
        }
    }

    // Flipping the vector (Used to have opposite order and signs)
    std::vector<double> reversed_z(z.rbegin(), z.rend());

    // Calculating LGL nodes based on z values
    for (int i = 0; i <= N; i++) {
        tnodes[i] = (reversed_z[i] + 1.0) / (2.0 / tf);
    }

    /////////////////////////////////////////////////////
    // Computing weights
    /////////////////////////////////////////////////////

    for (int i = 0; i <= N; i++) {
        w[i] = tf / (N * (N + 1) * std::pow(L[i][N], 2));
    }

    w[0] = tf / (N * (N + 1));
    w[N] = tf / (N * (N + 1));

    /////////////////////////////////////////////////////
    // Finding differentiation Matrix
    /////////////////////////////////////////////////////

    // Computing the differentiation matrix elements
    for (int i = 0; i <= N; i++) {
        for (int k = 0; k <= N; k++) {
            if (i == k) {
                Dm[i][k] = 0.0;
            } else {
                Dm[i][k] = (L[k][N] / L[i][N]) / (tnodes[k] - tnodes[i]);
            }

            if (i == 0 && k == 0) {
                Dm[i][k] = -N * (N + 1) / (2 * tf);
            }

            if (i == N && k == N) {
                Dm[i][k] = N * (N + 1) / (2 * tf);
            }
        }
    }
}

std::vector<double> LGL_PS::getNodes() {
    return tnodes;
}

std::vector<double> LGL_PS::getWeights() {
    return w;
}

std::vector<std::vector<double>> LGL_PS::getDifferentiationMatrix() {
    return Dm;
}
