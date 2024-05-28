#include "../../../../include/piecewiseintegrationmatrix.h"
#include <iostream>
#include <vector>

int main() {
    int K = 3;
    int N = 2;
    std::vector<double> T = {0, 1, 3, 5};
    
    std::vector<std::vector<double>> I = PiecewiseIntegrationMatrix(K, N, T);

    // Output the matrix for verification
    for (const auto& row : I) {
        for (double elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}


// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/test/ma_27/opt_hardcore_equi$ g++ -I~/dev/optimization/BeBOT_cpp_v2/include -o opt_hardcore_equi opt_hardcore_equi.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewiseintegrationmatrix.cpp
// ./opt_hardcore_equi 


