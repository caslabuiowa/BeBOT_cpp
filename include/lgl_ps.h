#ifndef LGL_PS_H
#define LGL_PS_H

#include <vector>

class LGL_PS {
public:
    LGL_PS(int N, double tf);
    void calculate();
    std::vector<double> getNodes();
    std::vector<double> getWeights();
    std::vector<std::vector<double>> getDifferentiationMatrix();

private:
    int N;
    double tf;
    std::vector<double> tnodes;
    std::vector<double> w;
    std::vector<std::vector<double>> Dm;
};

#endif
