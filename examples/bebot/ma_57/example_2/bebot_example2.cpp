#include "../../../../../ipopt/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../../ipopt/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/bebot.h"
#include "../../../../include/bernsteinpoly.h"
#include <iomanip>
#include "mkl.h"

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, double amax, double amin, double x10, double x1f, double x20, double x2f)
        : N_(N), tf_(tf), amax_(amax), amin_(amin), x10_(x10), x1f_(x1f), x20_(x20), x2f_(x2f), bebot_(N, tf_) {
        bebot_.calculate();
    }

    // write data to CSV
    void writeToCSV(const std::vector<double>& times, const std::vector<double>& values, const std::string& filename) {
        std::ofstream outFile(filename);

        if (!outFile.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        outFile << "Time,Value\n"; // CSV headers
        for (size_t i = 0; i < times.size(); ++i) {
            outFile << std::fixed << std::setprecision(6) << times[i] << "," << values[i] << "\n";
        }

        outFile.close();
    }

    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {
        // Define the number of variables, constraints, and Jacobian/Hessian non-zero elements.
        n = N_ + 2; 
        m = 2 * (N_ + 1);
        nnz_jac_g = n*m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Define variable and constraint bounds
        for (int i = 0; i < N_+ 2; i++) {
            x_l[i] = -std::numeric_limits<double>::infinity();
            x_u[i] = std::numeric_limits<double>::infinity();
            //std::cout << "x_l[" << i << "] = " << x_l[i] << std::endl;
            //std::cout << "x_u[" << i << "] = " << x_u[i] << std::endl;
        }

        x_l[0] = x10_;
        x_u[0] = x10_;
        x_l[(N_ + 2) - 2] = x1f_;
        x_u[(N_ + 2) - 2] = x1f_;
        x_l[(N_ + 2) - 1] = 0;
        x_u[(N_ + 2) - 1] = std::numeric_limits<double>::infinity();//std::numeric_limits<double>::infinity();;

        //std::cout << "x_l[" << 0 << "] = " << x_l[0] << std::endl;
        //std::cout << "x_u[" << 0 << "] = " << x_u[0] << std::endl;
        //std::cout << "x_l[" << (N_ + 2) - 2 << "] = " << x_l[(N_ + 2) - 2] << std::endl;
        //std::cout << "x_u[" << (N_ + 2) - 2 << "] = " << x_u[(N_ + 2) - 2] << std::endl;
        //std::cout << "x_l[" << (N_ + 2) - 1 << "] = " << x_l[(N_ + 2) - 1] << std::endl;
        //std::cout << "x_u[" << (N_ + 2) - 1 << "] = " << x_u[(N_ + 2) - 1] << std::endl;

        for (int i = 0; i < 2 * (N_ + 1); i++) {
            if (i == 0 || i == N_) {
                g_l[i] = 0;
                g_u[i] = 0;
                g_u[i] = std::numeric_limits<double>::infinity();
                //std::cout << "g_l[" << i << "] = " << g_l[i] << std::endl;
                //std::cout << "g_u[" << i << "] = " << g_u[i] << std::endl;
            } else {
                g_l[i] = amin_;
                g_u[i] = amax_;
                //std::cout << "g_l[" << i << "] = " << g_l[i] << std::endl;
                //std::cout << "g_u[" << i << "] = " << g_u[i] << std::endl;
            }
        }

        // Print constraint bounds
        //for (Index i = 0; i < 2 * (N_ + 1); ++i) {
        //    std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        //}

        return true;
    }

    // initialization of the starting point
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (int i = 0; i < N_ + 2; i++) {
            x[0] = 1;
            x[i] = 1;
            x[(N_+ 2) - 1] = 10;
            //std::cout << "x[" << i << "] = " << x[i] << std::endl;
            /*x[0] = -3;
            x[1] = -1.71;
            x[2] = 1.82;
            x[3] = 0.95;
            x[4] = -0.10;
            x[5] = 0;
            x[6] = 10;*/
        }
        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        // Objective function
        obj_value = x[(N_+ 2) - 1];
        //std::cout << "objective value = " << obj_value << std::endl;
        return true;
    }
    ///*
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Transpose the x vector to create a vertical column vector
        //std::vector<Number> x_transposed(N_ + 2);
        //for (Index i = 0; i < N_ + 1; i++) {
        //    x_transposed[i] = x[i];
        //    //std::cout << "x_transposed[" << i << "] = " << x_transposed[i] << std::endl;
        //}

        // Get the differentiation matrix
        Bebot Bebot(N_, x[(N_+ 2) - 1]);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();
        
        // Printing Dm matrix
        //std::cout << "Dm matrix:\n";
        //for (const auto& row : Dm) {
        //    for (const auto& element : row) {
        //        std::cout << element << " ";
        //    }
        //    std::cout << std::endl; 
        //}
        
        //std::cout << "N_ = " << N_ << std::endl;
        //std::cout << "x[" << (N_+ 2) - 1 << "] = " << x[(N_+ 2) - 1] << std::endl;

        // Perform the first operation g1 = (x' * Dm)'
        std::vector<Number> g1(N_ + 1, 0.0);
        //for (Index i = 0; i <= N_; i++) {
        //    for (Index j = 0; j <= N_; j++) {
                //std::cout << "_____start " << std::endl;
                //std::cout << "Dm " << Dm[i][j] << std::endl;
                //std::cout << "_____end " << Dm[i][j] << std::endl;
                //std::cout << "x[" << j << "] = " << x[j] << std::endl;
        //        g1[i] += x_transposed[j] * Dm[j][i];
        //        std::cout << "g1[" << i << "] = " << g[i] << std::endl;
         //   }
        //}

        ///*
        // preparing data for MKL in column-major order
        int N_mkl = N_ + 1; 
        //double* Dm_flat = new double[N_mkl * N_mkl]; 
        //double* x_flat = new double[N_mkl]; 
        double* result_flat = new double[N_mkl]; 

        // flattening the Dm matrix to Dm_flat in column-major order
        //for (int i = 0; i < N_mkl; ++i) {
        //    for (int j = 0; j < N_mkl; ++j) {
        //        Dm_flat[j * N_mkl + i] = Dm[i][j];
        //    }
        //}

        
        //for (int i = 0; i < N_mkl; ++i) {
        //    x_flat[i] = x[i];
        //}

        // After flattening Dm matrix to Dm_flat in column-major order
        //std::cout << "Dm_flat:" << std::endl;
        //for (int i = 0; i < N_mkl * N_mkl; ++i) {
        //    if (i % N_mkl == 0 && i != 0) std::cout << std::endl; // New line for each row
        //    std::cout << Dm_flat[i] << " ";
        //}
        //std::cout << std::endl << std::endl; // New line after printing the whole matrix

        // After copying x_transposed vector to x_flat
        //std::cout << "x_flat:" << std::endl;
        //for (int i = 0; i < N_mkl; ++i) {
        //    std::cout << x_flat[i] << " ";
        //}
        //std::cout << std::endl << std::endl; // New line after printing the vector


       
        cblas_dgemv(CblasColMajor, CblasTrans, N_mkl, N_mkl, 1.0, Dm.data(), N_mkl, x, 1, 0.0, result_flat, 1);

        // Print the result_flat array
        //std::cout << "result_flat after MKL multiplication:" << std::endl;
        //for (int i = 0; i < N_mkl; ++i) {
        //    std::cout << result_flat[i] << " ";
        //}
        //std::cout << std::endl << std::endl; 

        // Copy the result back into the std::vector or directly into the target structure
        for (int i = 0; i < N_mkl; ++i) {
            g[i] = result_flat[i]; 
            //std::cout << "g1_mkl[" << i << "] = " << g1[i] << std::endl;
        }

        // Use g1 for further operations as needed...

        // Cleanup
        //delete[] Dm_flat;
        //delete[] x_flat;
        //delete[] result_flat;

        //return true;
        //*/
        // Transpose g1 to create a vertical column vector
        //for (Index i = 0; i <= N_; i++) {
        //    g[i] = g1[i];
            //std::cout << "g[" << i << "] = " << g[i] << std::endl;
            
        //}

        // G2

        // Perform the second operation g2 = (g1' * Dm)'
        std::vector<Number> g2(N_ + 1, 0.0);
        //for (Index i = 0; i <= N_; i++) {
        //    for (Index j = 0; j <= N_; j++) {
        //        g2[i] += g1[j] * Dm[j][i];
        //    }
        //}


        // Prepare your data for MKL in column-major order again
        //double* g1_flat = new double[N_mkl]; // Allocate flat array for g1
        double* result_g2_flat = new double[N_mkl]; // Allocate flat array for the result of multiplying Dm by g1

        // Copy g1 vector to g1_flat
        //for (int i = 0; i < N_mkl; ++i) {
        //    g1_flat[i] = g1[i];
        //}

        // Perform the matrix-vector multiplication using MKL
        cblas_dgemv(CblasColMajor, CblasTrans, N_mkl, N_mkl, 1.0, Dm.data(), N_mkl, result_flat, 1, 0.0, result_g2_flat, 1);

        // Copy the result back into g2 vector or directly into the target structure
        //for (int i = 0; i < N_mkl; ++i) {
        //    g2[i] = result_g2_flat[i];
            //std::cout << "g2_mkl[" << i << "] = " << g2[i] << std::endl;
            // Optionally, print g2 values to verify
            // std::cout << "g2[" << i << "] = " << g2[i] << std::endl;
        //}

        // Cleanup
        //delete[] g1_flat;
        //delete[] result_g2_flat;

        // Transpose g2 to create the final vertical column vector
        for (Index i = 0; i <= N_; i++) {
            g[N_ + 1 + i] = result_g2_flat[i];
            //std::cout << "g[" << N_ + 1 + i << "] = " << g[N_ + 1 + i] << std::endl;
        }

        return true;
    }
    
    // Define the Jacobian of /the constraints
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values) {
        if (values == NULL) {
            // Return the structure of the Jacobian by setting iRow and jCol
            for (Index i = 0; i < m; i++) {
                for (Index j = 0; j < n; j++) {
                    iRow[i * n + j] = i;
                    jCol[i * n + j] = j;
                }
            }
        }
        return true;
    }

    // Define the gradient of the objective function
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
        return true;
    }

    // Method to finalize the solution
    virtual void finalize_solution(
        SolverReturn status, 
        Index n,
        const Number* x,
        const Number* z_L,
        const Number* z_U,
        Index m,
        const Number* g,
        const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq

    ) { solution_x_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x_[i] = x[i];
            //std::cout << "solution_x_[" << i << "] = " << solution_x_[i] << std::endl;
        }
        final_obj_value_ = obj_value;
        
        // Updating tf_ with the optimized value of tf (which is x[(N_+ 2) - 1])
        tf_ = x[(N_+ 2) - 1];
        
        // Now recalculating the BeBOT points with the updated final time tf_
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();
        
        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * obj_value / 999.0;
            //std::cout << "final_time_[" << i << "] = " << final_time_[i] << std::endl;
        }
        
        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        std::vector<std::vector<double>> solution_x_2d(1, std::vector<double>(solution_x_.begin(), solution_x_.end()));
        bernsteinpoly_result_ = BernsteinPoly(solution_x_2d, final_time_, 0, tf_);
        
        // After calculating final_time_ and bernsteinpoly_result_
        // Flatten bernsteinpoly_result_
        std::vector<double> flattened_result;
        for (const auto& row : bernsteinpoly_result_) {
            flattened_result.insert(flattened_result.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result, "x.csv");
        writeToCSV(bebot_.getNodes(), solution_x_, "x_controlpoints.csv");
        solution_x2_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x2_[i] = g[i];
            //std::cout << "solution_x2_[" << i << "] = " << solution_x2_[i] << std::endl;
        }
        std::vector<std::vector<double>> solution_x2_2d(1, std::vector<double>(solution_x2_.begin(), solution_x2_.end()));      
        bernsteinpoly_resultx2_ = BernsteinPoly(solution_x2_2d, final_time_, 0, tf_);
        
        std::vector<double> flattened_result1;
        for (const auto& row : bernsteinpoly_resultx2_) {
            flattened_result1.insert(flattened_result1.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result1, "x1.csv");
        writeToCSV(bebot_.getNodes(), solution_x2_, "x1_controlpoints.csv");
        solution_u_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_u_[i] = g[N_ + 1 + i];
            //std::cout << "solution_u_[" << i << "] = " << solution_u_[i] << std::endl;
        }
        std::vector<std::vector<double>> solution_u_2d(1, std::vector<double>(solution_u_.begin(), solution_u_.end()));        
        bernsteinpoly_resultu_ = BernsteinPoly(solution_u_2d, final_time_, 0, tf_);
        
        //std::cout << "solution_u_2d :" << std::endl;
        //for (const auto& row : bernsteinpoly_resultu_) {
        //    for (const auto& elem : row) {
        //        std::cout << elem << " ";
        //    }
        //    std::cout << std::endl; 
        //}

        std::vector<double> flattened_result2;
        for (const auto& row : bernsteinpoly_resultu_) {
            flattened_result2.insert(flattened_result2.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result2, "u.csv");
        writeToCSV(bebot_.getNodes(), solution_u_, "u_controlpoints.csv");
    }
    // Getter for the solution
    const std::vector<Number>& get_solution_x() const { return solution_x_; }

    // Getter for the final objective value
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    double tf_;
    double amax_;
    double amin_;
    double x10_;
    double x1f_;
    double x20_;
    double x2f_;
    Bebot bebot_;
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_;
    Number final_obj_value_;
    std::vector<double> final_time_;
    std::vector<std::vector<double>> bernsteinpoly_resultu_;
    std::vector<std::vector<double>> bernsteinpoly_resultx2_;
    std::vector<std::vector<double>> bernsteinpoly_result_;

public:
    const std::vector<std::vector<double>>& get_bernsteinpoly_result() const { 
        return bernsteinpoly_result_; }
};

int main() {
    int N = 5;
    double tf = 10.0;
    double amax = 5.0;
    double amin = -5.0;
    double x10 = -3.0;
    double x1f = 0.0;
    double x20 = 0.0;
    double x2f = 0.0;

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, tf, amax, amin, x10, x1f, x20, x2f);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetStringValue("linear_solver", "ma57");
    // A smaller number pivots for sparsity, a larger number pivots for stability
    //app->Options()->SetNumericValue("ma57_pivtol", 1e-8);//0.99 // 1e-8 // between 0 and 1
    // Ipopt may increase pivtol as high as ma27_pivtolmax to get a more accurate solution to the linear system
    //app->Options()->SetNumericValue("ma57_pivtolmax", 0.99);//0.99 // 0.0001 // between 0 and 1
    // The initial integer workspace memory = liw_init_factor * memory required by unfactored system. 
    // Ipopt will increase the workspace size by ma27_meminc_factor if required.
    //app->Options()->SetNumericValue("ma57_liw_init_factor", 5.0); // 5.0 has to be 
    // The initial real workspace memory = la_init_factor * memory required by unfactored system. 
    // Ipopt will increase the workspace size by ma27_meminc_factor if required
    //app->Options()->SetNumericValue("ma57_la_init_factor", 5.0); // 5.0
    // If the integer or real workspace is not large enough, Ipopt will increase its size by this factor.
    //app->Options()->SetNumericValue("ma57_meminc_factor", 5.0); // 5.0

    app->Options()->SetStringValue("mu_strategy", "adaptive");
    
    app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
    app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    // Set the Hessian approximation method to limited-memory
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // Adjust the maximum number of iterations
    app->Options()->SetIntegerValue("max_iter", 5000); // Change to my desired maximum iterations

    // Adjust the convergence tolerance
    app->Options()->SetNumericValue("tol", 1e-6); // Change to my desired tolerance



    app->RethrowNonIpoptException(true);
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << "IPOPT initialization failed!" << std::endl;
        return -1;
    }

    status = app->OptimizeTNLP(pointSetProblem);

    if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
        // Retrieve the optimal solution and objective value from the problem
        const auto& solution_x = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_solution_x();
        Number final_obj_value = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_final_obj_value();
    
        std::cout << "Optimal Solution (x): ";
        for (Index i = 0; i < solution_x.size(); i++) {
            std::cout << solution_x[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Optimal Objective Value: " << final_obj_value << std::endl;
        
        // Retrieve the Bernstein polynomial result from the problem
        const auto& bernsteinpoly_result = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_bernsteinpoly_result();
        //std::cout << "Bernstein Polynomial Result (xN): " << std::endl;
        //for (const auto& row : bernsteinpoly_result) {
        //    for (double value : row) {
        //        std::cout << value << " ";
        //    }
        //    std::cout << std::endl; // New line for each row
        //}
    } else {
        std::cout << "IPOPT optimization failed with status " << status << std::endl;
    }

    return 0;
}

// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_2$ g++ -o bebot_example2_ma57_mkl ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_2/bebot_example2.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_2$ export LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_2$ ./bebot_example2_ma57_mkl