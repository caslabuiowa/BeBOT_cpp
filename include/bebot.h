#ifndef BEBOT_H
#define BEBOT_H

#include <vector>

class Bebot {
public:
    Bebot(int N, double T);
    void calculate();
    std::vector<double> getNodes();
    std::vector<double> getWeights();
    const std::vector<double>& getDifferentiationMatrix() const;

private:
    int N;
    double T;
    std::vector<double> tnodes;
    std::vector<double> w;
    std::vector<double> Dm;
};

#endif
