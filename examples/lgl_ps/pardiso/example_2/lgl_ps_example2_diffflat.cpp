#include "../../../../Ipopt_pardiso/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_pardiso/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/lgl_ps.h"
#include "../../../../include/lagrangepoly.h"
#include <iomanip>

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, double amax, double amin, double x10, double x1f, double x20, double x2f)
        : N_(N), tf_(tf), amax_(amax), amin_(amin), x10_(x10), x1f_(x1f), x20_(x20), x2f_(x2f), lgl(N, tf_) {
        lgl.calculate();
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
                //std::cout << "g[" << i << "] = " << g[i] << std::endl;
            } else if (i >= 1 && i < N_ + 1) {
                g_l[i] = -std::numeric_limits<double>::infinity();
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
        std::vector<Number> x_transposed(N_ + 2);
        for (Index i = 0; i < N_ + 1; i++) {
            x_transposed[i] = x[i];
            //std::cout << "x_transposed[" << i << "] = " << x_transposed[i] << std::endl;
        }

        // Get the differentiation matrix
        LGL_PS lgl(N_, x[(N_+ 2) - 1]);
        lgl.calculate();
        std::vector<std::vector<double>> Dm = lgl.getDifferentiationMatrix();
        
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
        for (Index i = 0; i <= N_; i++) {
            for (Index j = 0; j <= N_; j++) {
                //std::cout << "_____start " << std::endl;
                //std::cout << "Dm " << Dm[i][j] << std::endl;
                //std::cout << "_____end " << Dm[i][j] << std::endl;
                //std::cout << "x[" << j << "] = " << x[j] << std::endl;
                g1[i] += x_transposed[j] * Dm[j][i];
            }
        }

        // Transpose g1 to create a vertical column vector
        for (Index i = 0; i <= N_; i++) {
            g[i] = g1[i];
            //std::cout << "g[" << i << "] = " << g[i] << std::endl;
            
        }

        // Perform the second operation g2 = (g1' * Dm)'
        std::vector<Number> g2(N_ + 1, 0.0);
        for (Index i = 0; i <= N_; i++) {
            for (Index j = 0; j <= N_; j++) {
                g2[i] += g1[j] * Dm[j][i];
            }
        }

        // Transpose g2 to create the final vertical column vector
        for (Index i = 0; i <= N_; i++) {
            g[N_ + 1 + i] = g2[i];
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
        
        // Now recalculate the LGL points with the updated final time tf_
        lgl = LGL_PS(N_, tf_); // You might need to modify LGL_PS to allow updating N and tf
        lgl.calculate();

        // Calculating final time t using obj_value as final tf for lagrangepoly library
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * obj_value / 999.0;
            //std::cout << "final_time_[" << i << "] = " << final_time_[i] << std::endl;
        }
        lagrange_result_ = LagrangePoly(solution_x_, lgl.getNodes(), final_time_);
        //std::cout << "Finalizing solution, preparing to write CSV file." << std::endl;

        // After calculating final_time_ and lagrange_result_
        writeToCSV(final_time_, lagrange_result_, "x.csv");
        writeToCSV(lgl.getNodes(), solution_x_, "x_controlpoints.csv");

        //std::cout << "CSV file writing completed." << std::endl;
        
        // Print tnodes here
        //const std::vector<double>& tnodes = lgl.getNodes();
        //std::cout << "tnodes: ";
        //for (double tn : tnodes) {
        //    std::cout << tn << " ";
        //}
        //std::cout << std::endl;

        solution_x2_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x2_[i] = g[i];
            //std::cout << "solution_x2_[" << i << "] = " << solution_x2_[i] << std::endl;
        }

        lagrange_resultx2_ = LagrangePoly(solution_x2_, lgl.getNodes(), final_time_);
        //std::cout << "Finalizing solution, preparing to write CSV file." << std::endl;
        // After calculating final_time_ and lagrange_result_
        writeToCSV(final_time_, lagrange_resultx2_, "x1.csv");
        writeToCSV(lgl.getNodes(), solution_x2_, "x1_controlpoints.csv");
        //std::cout << "CSV file writing completed." << std::endl;

        solution_u_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_u_[i] = g[N_ + 1 + i];
            //std::cout << "solution_u_[" << i << "] = " << solution_u_[i] << std::endl;
        }

        lagrange_resultu_ = LagrangePoly(solution_u_, lgl.getNodes(), final_time_);
        //std::cout << "Finalizing solution, preparing to write CSV file." << std::endl;
        // After calculating final_time_ and lagrange_result_
        writeToCSV(final_time_, lagrange_resultu_, "u.csv");
        writeToCSV(lgl.getNodes(), solution_u_, "u_controlpoints.csv");
        //std::cout << "CSV file writing completed." << std::endl;
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
    LGL_PS lgl;
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_; 
    Number final_obj_value_; 
    std::vector<double> final_time_;
    std::vector<double> lagrange_resultu_;
    std::vector<double> lagrange_resultx2_;
    std::vector<double> lagrange_result_;

public:
    const std::vector<double>& get_lagrange_result() const { return lagrange_result_; }

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

    //LGL_PS lgl(N, tf);
    //lgl.calculate();

    //std::vector<double> tnodes = lgl.getNodes();
    //std::vector<double> w = lgl.getWeights();
    //std::vector<std::vector<double>> Dm = lgl.getDifferentiationMatrix();

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, tf, amax, amin, x10, x1f, x20, x2f);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetStringValue("linear_solver", "pardiso");
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
    app->Options()->SetNumericValue("tol", 1e-2); // Change to my desired tolerance



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
        
        // Retrieve the Lagrange polynomial result from the problem
        const auto& lagrange_result = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_lagrange_result();
        //std::cout << "Lagrange Polynomial Result (xN): ";
        for (double value : lagrange_result) {
            //std::cout << value << " ";
        }
        std::cout << std::endl; // Make sure to end the line after printing
    } else {
        std::cout << "IPOPT optimization failed with status " << status << std::endl;
    }

    return 0;
}

// vpetrov@lnx-me002:/local/vol00/home/vpetrov/dev/optimization/BeBOT_cpp_v2/examples/lgl_ps/pardiso/example_2$ g++ -o lgl_ps_example2_diffflat_v2 ~/dev/optimization/BeBOT_cpp_v2/examples/lgl_ps/pardiso/example_2/lgl_ps_example2_diffflat.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/lgl_ps.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/lagrangepoly.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I/usr/local/include -L/usr/local/lib -lipopt -L/usr/local/lib/pardiso/panua-pardiso-20230908-linux/lib -lpardiso -ldl -lm -lstdc++
// vpetrov@lnx-me002:/local/vol00/home/vpetrov/dev/optimization/BeBOT_cpp_v2/examples/lgl_ps/pardiso/example_2$ export LD_LIBRARY_PATH=/usr/local/lib/pardiso/panua-pardiso-20230908-linux/lib:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:/local/vol00/home/vpetrov/dev/optimization/BeBOT_cpp_v2/examples/lgl_ps/pardiso/example_2$ ./lgl_ps_example2_diffflat_v2
