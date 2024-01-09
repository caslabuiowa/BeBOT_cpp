#ifndef PIECEWISEBEBOT_H
#define PIECEWISEBEBOT_H

#include <vector>

class PiecewiseBeBOT {
public:
    PiecewiseBeBOT(int N, const std::vector<double>& tknots);
    void calculate();
    std::vector<double> getNodes() const;
    std::vector<double> getWeights() const;
    std::vector<std::vector<double>> getDifferentiationMatrix() const;

private:
    int N;
    std::vector<std::vector<double>> transformedTknots;
    std::vector<double> tnodes;
    std::vector<double> w;
    std::vector<std::vector<double>> Dm;
    std::vector<double> originalTknots;

    void transformTknots(); // New method for transforming tknots
};

#endif
