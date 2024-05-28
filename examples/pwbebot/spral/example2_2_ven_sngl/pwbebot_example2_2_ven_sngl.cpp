#include "../../../../Ipopt_spral_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_spral_solver/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/piecewisebebot.h"
#include "../../../../include/piecewisebernsteinpoly.h"
#include <iomanip>

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, int M, double tf, double umax, double umin, double y0, double yf)
        : N_(N), M_(M), tf_(tf), umax_(umax), umin_(umin), y0_(y0), yf_(yf),
          piecewiseBebot_(N, generateTknots()) { // Initialize PiecewiseBeBOT in the initializer list
        piecewiseBebot_.calculate();
    }

    // Write data to CSV
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
        n = M_* (N_ + 1); 
        m = M_ * (N_ + 1) + (M_ + 1);
        nnz_jac_g = n*m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Define variable bounds
        for (Index i = 0; i < M_ * (N_ + 1); i++) {
            x_l[i] = -std::numeric_limits<double>::infinity();
            x_u[i] = std::numeric_limits<double>::infinity();
        }

        // First element bounds
        x_l[0] = y0_;
        x_u[0] = y0_;

        // Continuity constraints between segments for x
        for (Index segment = 0; segment < M_; ++segment) {
            Index endOfSegment = (segment + 1) * (N_ + 1) - 1;
            Index startOfNextSegment = (segment + 1) * (N_ + 1);

            x_l[endOfSegment] = x_l[startOfNextSegment];
            x_u[endOfSegment] = x_u[startOfNextSegment];
        }

        // Last element bounds
        x_l[M_ * (N_ + 1) - 1] = yf_;
        x_u[M_ * (N_ + 1) - 1] = yf_;

        // Time variable bounds
        //x_l[M_ * (N_ + 1)] = 0;
        //x_u[M_ * (N_ + 1)] = std::numeric_limits<double>::infinity();


        // Print the bounds for each variable
        for (Index i = 0; i < M_ * (N_ + 1); ++i) {
            std::cout << "x_l[" << i << "] = " << x_l[i] << ", x_u[" << i << "] = " << x_u[i] << std::endl;
        }

        // Define constraint bounds
        for (int i = 0; i < M_ * (N_ + 1) + (M_ + 1); i++) {
            if (i >= M_ * (N_ + 1)) {
                // For indices after M_ * (N_ + 1), set g_l and g_u between 0 and tf_
                g_l[i] = 0;
                g_u[i] = tf_;
            } else {
                // For indices before M_ * (N_ + 1), set g_l and g_u between -umin_ and umax_
                g_l[i] = umin_;
                g_u[i] = umax_;
            }

            g_l[M_ * (N_ + 1) + (M_ + 1)-1] = 2;
            g_u[M_ * (N_ + 1) + (M_ + 1)-1] = 2;
   
        }

        // Print constraint bounds
        for (Index i = 0; i < M_ * (N_ + 1) + (M_+ 1); ++i) {
            std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        }

        return true;
    }

    // initialization of the starting point
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (int i = 0; i < M_ * (N_ + 1); i++) {
            x[i] = 1;
        }

        //for (Index i = 0; i < M_ * (N_ + 1) + 1; ++i) {
        //    std::cout << "x_in_guess[" << i << "] = " << x[i] << std::endl;
        //}

        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        // Objective function
        //obj_value = x[M_ * (N_ + 1)];
        if (g_values_.empty() || g_values_.size() != M_ * (N_ + 1) + (M_ + 1)) {
            std::cerr << "Error: g_values_ is not correctly populated." << std::endl;
            return false;
        }

        std::vector<double> weights = piecewiseBebot_.getWeights();

        obj_value = 0;
        for (Index i = 0; i < n; ++i) {
            obj_value += weights[i] * (3 * g_values_[i] - 2 * x[i]);
        }

        //std::cout << "objective value = " << obj_value << std::endl;
        return true;
    }
    ///*
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Transpose the x vector to create a vertical column vector
        std::vector<Number> x_transposed(M_ * (N_ + 1));
        for (Index i = 0; i < M_ * (N_ + 1); i++) {
            x_transposed[i] = x[i];
            std::cout << "x[" << i << "] = " << x[i] << std::endl;
            //std::cout << "x_transposed[" << i << "] = " << x_transposed[i] << std::endl;
        }
        
        // Get the differentiation matrix
        //tf_ = x[M_ * (N_ + 1)];
        std::vector<double> tknots = generateTknots();

        // Print tknots
        //std::cout << "tknots: ";
        //for (const auto& knot : tknots) {
        //    std::cout << knot << " ";
        //}
        //std::cout << std::endl;

        
        PiecewiseBeBOT piecewiseBebot(N_, tknots);
        piecewiseBebot.calculate();
        std::vector<std::vector<double>> Dm = piecewiseBebot.getDifferentiationMatrix();
        
        // Printing Dm matrix
        //std::cout << "Dm matrix:\n";
        //for (const auto& row : Dm) {
        //    for (const auto& element : row) {
        //        std::cout << element << " ";
        //    }
        //    std::cout << std::endl; 
        //}
        
        //std::cout << "N_ = " << N_ << std::endl;
        //std::cout << "x[" << M_ * (N_ + 1) << "] = " << x[M_ * (N_ + 1)] << std::endl;

        // Perform the first operation g1 = (x' * Dm)'
        std::vector<Number> g1(M_ * (N_ + 1), 0.0);  // Resizing g1 to handle all segments

        for (Index segment = 0; segment < M_; ++segment) {
            for (Index i = 0; i <= N_; i++) {
                for (Index j = 0; j <= N_; j++) {
                    // Computing for each segment separately
                    Index xIndex = segment * (N_ + 1) + j;
                    Index g1Index = segment * (N_ + 1) + i;
                    //std::cout << "x_transposed[" << xIndex << "] = " << x_transposed[xIndex] << std::endl;

                    g1[g1Index] += x_transposed[xIndex] * Dm[j][i];
                }
            }
        }


        // Print g1 values after calculation
        //std::cout << "g1 vector after calculation:" << std::endl;
        //for (Index i = 0; i < g1.size(); i++) {
        //    std::cout << "g1_afc[" << i << "] = " << g1[i] << std::endl;
        //}
    

        // Continuity constraints for g1
        //for (Index segment = 0; segment < M_ - 1; ++segment) {
        //    Index endOfSegment = (segment + 1) * (N_ + 1) - 1;
        //    Index startOfNextSegment = (segment + 1) * (N_ + 1);
        //    g1[endOfSegment] = g1[startOfNextSegment];
        //}

        // Print g1 values after applying continuity constraints
        //std::cout << "g1 vector after applying continuity constraints:" << std::endl;
        //for (Index i = 0; i < g1.size(); i++) {
        //    std::cout << "g1_afc_c[" << i << "] = " << g1[i] << std::endl;
        //}

        // Print g1 values
        //for (Index i = 0; i < g1.size(); ++i) {
        //    std::cout << "g1[" << i << "] = " << g1[i] << std::endl;
        //}

        // Transpose g1 to create a vertical column vector
        for (Index i = 0; i < g1.size(); i++) {
            g[i] = g1[i];
            //std::cout << "g[" << i << "] = " << g[i] << std::endl;
            
        }

        
        Index offset = M_ * (N_ + 1);

        if (M_ > 1) {
            for (Index i = 0; i < M_ - 1; ++i) {
                Index endOfSegment = (i + 1) * (N_ + 1) - 1;
                Index startOfNextSegment = (i + 1) * (N_ + 1);

                g[offset + i] = x[endOfSegment] - x[startOfNextSegment]; // Continuity constraint for x
                
                Index g_continuity_index = offset + M_ - 1 + i;
                g[g_continuity_index] = g1[endOfSegment] - g1[startOfNextSegment]; // Continuity constraint for g

            }
        }
    

        // Define additional constraints for g values after M_ * (N_ + 1)
        Number increment = tf_ / M_; // Increment value to ensure constraints
        for (Index i = offset + M_ - 1; i < offset + M_ + 1; ++i) {
            g[i] = (i - offset - M_ + 2) * increment; // Incrementally increasing each g element
        }

        // Ensure all values are within 0 to 2
        for (Index i = offset + M_ - 1; i < offset + M_ + 1; ++i) {
            if (g[i] < 0) g[i] = 0;
            if (g[i] > 2) g[i] = 2;
        }


        // Transpose g2 to create the final vertical column vector
        //for (Index i = 0; i < g2.size(); i++) {
        //    g[N_ + 1 + i] = g2[i];
        //    std::cout << "g[" << N_ + 1 + i << "] = " << g[N_ + 1 + i] << std::endl;
        //}

        
        
        //std::cout << "Complete g vector:" << std::endl;
        for (Index i = 0; i < M_ * (N_ + 1) + (M_ + 1); ++i) {
            std::cout << "g[" << i << "] = " << g[i] << std::endl;
        }

        // Store g values
        g_values_.assign(g, g + m);

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

    ) 
    
    
    { 

    }
    // Getter for the solution
    const std::vector<Number>& get_solution_x() const { return solution_x_; }

    // Getter for the final objective value
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    int M_;
    double tf_;
    double umax_;
    double umin_;
    double y0_;
    double yf_;

    PiecewiseBeBOT piecewiseBebot_; 
    std::vector<Number> g_values_;
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_;
    Number final_obj_value_;
    std::vector<double> final_time_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_resultu_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_resultx2_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_;

    // Helper function to generate tknots
    std::vector<double> generateTknots() {
        std::vector<double> tknots;
        //std::cout << "tf_ = " << tf_ << std::endl;
        //std::cout << "M = " << M_ << std::endl;
        double interval = tf_ / M_;
        for (int i = 0; i <= M_; ++i) {
            tknots.push_back(i * interval);
        }
        return tknots;
    }

public:
    
};

int main() {
    int N = 5;
    int M = 2;
    double tf = 2;
    double umax = 2.0;
    double umin = -2.0;
    double y0 = 4.0;
    double yf = 39.392;

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, M, tf, umax, umin, y0, yf);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    //app->Options()->SetStringValue("linear_solver", "custom");
    app->Options()->SetStringValue("linear_solver", "spral");
    //app->Options()->SetStringValue("custom_linear_solver_library", "libspral.so");
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
    app->Options()->SetIntegerValue("max_iter", 10); // Change to my desired maximum iterations

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
        // Process optimization results here
        const auto& solution_x = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_solution_x();
        Number final_obj_value = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_final_obj_value();
    
        std::cout << "Optimal Solution (x): ";
        for (Index i = 0; i < solution_x.size(); i++) {
            std::cout << solution_x[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Optimal Objective Value: " << final_obj_value << std::endl;
    } else {
        std::cout << "IPOPT optimization failed with status " << status << std::endl;
    }

    return 0;
}

// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example_2$ g++ -o pwbebot_example2_2_ven_sngl ~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example2_2_ven_sngl/pwbebot_example2_2_ven_sngl.cpp ~/dev/optimization/BeBOT_cpp/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -lspral -ldl -lm -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example_2$ export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example_2$ export OMP_CANCELLATION=TRUE
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example_2$ ./pwbebot_example2_2_ven_sngl


//vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example2_2_ven_sngl$ g++ -o pwbebot_example2_2_ven_sngl ~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example2_2_ven_sngl/pwbebot_example2_2_ven_sngl.cpp ~/dev/optimization/BeBOT_cpp/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -lspral -ldl -lm -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example_2$ export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example_2$ export OMP_CANCELLATION=TRUE
//vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp/examples/pwbebot/spral/example2_2_ven_sngl$ ./pwbebot_example2_2_ven_sngl 
