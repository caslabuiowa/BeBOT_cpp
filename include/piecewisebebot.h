#ifndef PIECEWISEBEBOT_H
#define PIECEWISEBEBOT_H

#include <vector>

class PiecewiseBeBOT {
public:
    PiecewiseBeBOT(int N, const std::vector<double>& tknots);
    void calculate();
    std::vector<double> getNodes() const;
    std::vector<double> getWeights() const;
    std::vector<std::vector<double>> getDifferentiationMatrix() const; // test, can delete bc it's no longer needed
    const std::vector<double>& getDifferentiationMatrixFlat() const; 

private:
    int N;
    std::vector<std::vector<double>> transformedTknots;
    std::vector<double> tnodes;
    std::vector<double> w;
    std::vector<double> Dm_flat; 
    std::vector<double> originalTknots;

    void transformTknots(); 
};

#endif 
