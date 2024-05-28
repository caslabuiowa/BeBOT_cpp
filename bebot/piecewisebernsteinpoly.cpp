#include "../include/piecewisebernsteinpoly.h"
#include "../include/bernsteinpoly.h"
#include <iostream>

std::vector<std::vector<double>> PiecewiseBernsteinPoly(const std::vector<std::vector<double>>& Cp, 
                                                        const std::vector<double>& tknots, 
                                                        const std::vector<double>& t) {
    size_t M = tknots.size() - 1;
    //std::cout << "M: " << M << std::endl;
    size_t dim = Cp.size();
    //std::cout << "dim: " << dim << std::endl;
    size_t totN = Cp[0].size();
    //std::cout << "totN: " << totN << std::endl;
    size_t N = totN / M - 1;
    //std::cout << "M: " << M << std::endl;

    std::vector<std::vector<double>> poly_t(dim, std::vector<double>(t.size(), 0.0));

    for (size_t i = 0; i < M; ++i) {
        size_t start_idx = i * (N + 1);
        //std::cout << "start_idx: " << start_idx << std::endl;
        size_t end_idx = (i < M - 1) ? start_idx + N + 1 : totN;
        //std::cout << "end_idx: " << end_idx << std::endl;

        for (size_t k = 0; k < t.size(); ++k) {
            bool isInInterval = (i < M - 1) ? (t[k] >= tknots[i] && t[k] < tknots[i + 1]) : (t[k] >= tknots[i] && t[k] <= tknots[i + 1]);
            //std::cout << "isInInterval: " << isInInterval << std::endl;
            if (isInInterval) {
                std::vector<std::vector<double>> Cp_segment(dim, std::vector<double>(N + 1));
                for (size_t d = 0; d < dim; ++d) {
                    std::copy(Cp[d].begin() + start_idx, Cp[d].begin() + end_idx, Cp_segment[d].begin());
                }

                std::vector<double> single_time = {t[k]};
                std::vector<std::vector<double>> result = BernsteinPoly(Cp_segment, single_time, tknots[i], tknots[i + 1]);
                for (size_t d = 0; d < dim; ++d) {
                    poly_t[d][k] = result[d][0];
                    //std::cout << "poly_t[" << d << "][" << k << "] = " << poly_t[d][k] << std::endl;

                }
            }
        }
    }

    return poly_t;
}
