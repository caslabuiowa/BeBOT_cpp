#include "../../../../include/piecewisebernsteinproduct.h"
#include <iostream>
#include <vector>

int main() {
    std::vector<double> A = {0.178056961080680, 0.145083011153051, 0.145083011153051, 0.112109061225423, 0.112109061225423, 0.0729303173154614, 0.0729303173154614, 0.0368520223781098, 0.0368520223781098, 9.68552570090875e-08, 9.68552570090875e-08, -0.0368528469536680, -0.0368528469536680, -0.0729290223440103, -0.0729290223440103, -0.112109600067686, -0.112109600067686, -0.145083160068060, -0.145083160068060, -0.178056720068434};
    std::vector<double> B = A;  
    int K = 10;
    int N = 1;

    std::vector<double> Cp = PiecewiseBernsteinProduct(A, B, K, N);
    for (double val : Cp) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}

// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/test/ma_27/opt_hardcore_equi$ g++ -I~/dev/optimization/BeBOT_cpp_v2/include -o testbp testbp.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebernsteinproduct.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinproduct.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp
