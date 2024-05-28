#include "../../../../Ipopt_pardiso/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_pardiso/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/piecewisebebot.h"
#include "../../../../include/piecewisebernsteinpoly.h"
#include <iomanip>
#include "mkl.h"


using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, int M, double tf, double amax, double amin, double x10, double x1f, double x20, double x2f)
        : N_(N), M_(M), tf_(tf), amax_(amax), amin_(amin), x10_(x10), x1f_(x1f), x20_(x20), x2f_(x2f),
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
        n = M_* (N_ + 1) + 1; 
        m = M_* 2 * (N_ + 1) + 2*(M_ - 1);
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
        x_l[0] = x10_;
        x_u[0] = x10_;

        // Continuity constraints between segments for x
        for (Index segment = 0; segment < M_; ++segment) {
            Index endOfSegment = (segment + 1) * (N_ + 1) - 1;
            Index startOfNextSegment = (segment + 1) * (N_ + 1);

            x_l[endOfSegment] = x_l[startOfNextSegment];
            x_u[endOfSegment] = x_u[startOfNextSegment];
        }

        // Last element bounds
        x_l[M_ * (N_ + 1) - 1] = x1f_;
        x_u[M_ * (N_ + 1) - 1] = x1f_;

        // Time variable bounds
        x_l[M_ * (N_ + 1)] = 0;
        x_u[M_ * (N_ + 1)] = std::numeric_limits<double>::infinity();


        // Print the bounds for each variable
        //for (Index i = 0; i < M_ * (N_ + 1) + 1; ++i) {
        //    std::cout << "x_l[" << i << "] = " << x_l[i] << ", x_u[" << i << "] = " << x_u[i] << std::endl;
        //}

        // Define constraint bounds
        for (int i = 0; i < M_ * 2 * (N_ + 1); i++) {
            g_l[i] = -std::numeric_limits<double>::infinity();
            g_u[i] = std::numeric_limits<double>::infinity();

            // Set bounds for the first and last element of the first half of g
            if (i == 0 || i == M_ * (N_ + 1) - 1) {
                g_l[i] = 0;
                g_u[i] = 0;
            }

            // Second half (M * 2 * (N + 1))
            if (i >= M_ * (N_ + 1)) {
                g_l[i] = amin_;
                g_u[i] = amax_;
            }

            Index offset = M_ * 2 * (N_ + 1); // Assuming g already has 2*M*(N+1) bounds set

            if (M_ > 1) {
                for (Index i = 0; i < 2* (M_ - 1); ++i) {
                    g_l[offset + i] = 0; // Lower bound for continuity constraint
                    g_u[offset + i] = 0; // Upper bound for continuity constraint
                }
            }
        }

        // Print constraint bounds
        //for (Index i = 0; i < M_ * 2 * (N_ + 1) + M_-1; ++i) {
        //    std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        //}

        return true;
    }

    // initialization of the starting point
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (int i = 0; i < M_ * (N_ + 1); i++) {
            x[0] = 1;
            x[i] = 1;
            x[M_ * (N_ + 1)] = 10;
        }

        //for (Index i = 0; i < M_ * (N_ + 1) + 1; ++i) {
        //    std::cout << "x[" << i << "] = " << x[i] << std::endl;
        //}

        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        // Objective function
        obj_value = x[M_ * (N_ + 1)];
        //std::cout << "objective value = " << obj_value << std::endl;
        return true;
    }
    ///*
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Transpose the x vector to create a vertical column vector
        //std::vector<Number> x_transposed(M_ * (N_ + 1));
        //for (Index i = 0; i <= M_ * (N_ + 1); i++) {
        //    x_transposed[i] = x[i];
            //std::cout << "x[" << i << "] = " << x[i] << std::endl;
            //std::cout << "x_transposed[" << i << "] = " << x_transposed[i] << std::endl;
        //}
        
        // Get the differentiation matrix
        tf_ = x[M_ * (N_ + 1)];
        std::vector<double> tknots = generateTknots();

        // Print tknots
        //std::cout << "tknots: ";
        //for (const auto& knot : tknots) {
        //    std::cout << knot << " ";
        //}
        //std::cout << std::endl;

        
        PiecewiseBeBOT piecewiseBebot(N_, tknots );
        piecewiseBebot.calculate();
        //std::vector<std::vector<double>> Dm = piecewiseBebot.getDifferentiationMatrix();
        
        // Printing Dm matrix
        //std::cout << "Dm matrix:\n";
        //for (const auto& row : Dm) {
        //    for (const auto& element : row) {
        //        std::cout << element << " ";
        //    }
        //    std::cout << std::endl; 
        //}
        const std::vector<double>& Dm_flat = piecewiseBebot.getDifferentiationMatrixFlat();
        //std::cout << "N_ = " << N_ << std::endl;
        //std::cout << "x[" << M_ * (N_ + 1) << "] = " << x[M_ * (N_ + 1)] << std::endl;

        // Perform the first operation g1 = (x' * Dm)'
        std::vector<Number> g1(M_ * (N_ + 1), 0.0);  // Resizing g1 to handle all segments
        for (int segment = 0; segment < M_; ++segment) {
            cblas_dgemv(CblasRowMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, &Dm_flat[segment * (N_ + 1) * (N_ + 1)],
                        N_ + 1, &x[segment * (N_ + 1)], 1, 0.0, &g1[segment * (N_ + 1)], 1);
        }


        //for (Index segment = 0; segment < M_; ++segment) {
        //    for (Index i = 0; i <= N_; i++) {
        //        for (Index j = 0; j <= N_; j++) {
        //            // Computing for each segment separately
        //            Index xIndex = segment * (N_ + 1) + j;
        //            Index g1Index = segment * (N_ + 1) + i;
        //            //std::cout << "x_transposed[" << xIndex << "] = " << x_transposed[xIndex] << std::endl;

        //            g1[g1Index] += x_transposed[xIndex] * Dm[j][i];
        //        }
        //    }
        //}


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

        // Perform the second operation g2 = (g1' * Dm)'
        //std::vector<Number> g2(M_ * (N_ + 1), 0.0);  // Resizing g2 to handle all segments
        //for (Index segment = 0; segment < M_; ++segment) {
        //    for (Index i = 0; i <= N_; i++) {
        //        for (Index j = 0; j <= N_; j++) {
        //            // Computing for each segment separately
        //            Index g1Index = segment * (N_ + 1) + j;
        //            Index g2Index = segment * (N_ + 1) + i;
        //            g2[g2Index] += g1[g1Index] * Dm[j][i];
        //        }
        //    }
        //}
        
        std::vector<Number> g2(M_ * (N_ + 1), 0.0);

        for (int segment = 0; segment < M_; ++segment) {
            cblas_dgemv(CblasRowMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, &Dm_flat[segment * (N_ + 1) * (N_ + 1)],
                        N_ + 1, &g1[segment * (N_ + 1)], 1, 0.0, &g2[segment * (N_ + 1)], 1);
        }
        
        
        Index offset = M_ * 2 * (N_ + 1);

        if (M_ > 1) {
        for (Index i = 0; i < M_ - 1; ++i) {
            Index endOfSegment = (i + 1) * (N_ + 1) - 1;
            Index startOfNextSegment = (i + 1) * (N_ + 1);

            g[offset + i] = x[endOfSegment] - x[startOfNextSegment]; // Continuity constraint for x
            
            Index g_continuity_index = offset + M_ - 1 + i;
            g[g_continuity_index] = g1[endOfSegment] - g1[startOfNextSegment]; // Continuity constraint for g

        }
    }
    


        // Transpose g2 to create the final vertical column vector
        //for (Index i = 0; i < g2.size(); i++) {
        //    g[N_ + 1 + i] = g2[i];
        //    std::cout << "g[" << N_ + 1 + i << "] = " << g[N_ + 1 + i] << std::endl;
        //}

        // Transpose g2 to create the final vertical column vector
        for (Index segment = 0; segment < M_; ++segment) {
            for (Index i = 0; i <= N_; i++) {
                Index gIndex = M_ * (N_ + 1) + segment * (N_ + 1) + i;  // Corrected index calculation
                g[gIndex] = g2[segment * (N_ + 1) + i];
                //std::cout << "g[" << gIndex << "] = " << g[gIndex] << std::endl;
            }
        }
        
        //std::cout << "Complete g vector:" << std::endl;
        //for (Index i = 0; i < M_ * 2 * (N_ + 1) +M_-1; ++i) {
        //    std::cout << "g[" << i << "] = " << g[i] << std::endl;
        //}

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
    
    // x1 - x
    { solution_x_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x_[i] = x[i];
            //std::cout << "solution_x_[" << i << "] = " << solution_x_[i] << std::endl;
        }
        final_obj_value_ = obj_value;
        
        // Updating tf_ with the optimized value of tf (which is x[M_ * (N_ + 1))
        tf_ = x[M_ * (N_ + 1)];
        //std::cout << "tf_ = " << tf_ << std::endl;

        // Recalculate the PiecewiseBeBOT points with the updated final time tf_
        std::vector<double> tknots = generateTknots();
        piecewiseBebot_ = PiecewiseBeBOT(N_, tknots);
        piecewiseBebot_.calculate();

    
        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * obj_value / 999.0;
            //std::cout << "final_time_[" << i << "] = " << final_time_[i] << std::endl;
        }

        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        std::vector<std::vector<double>> solution_x_2d(1, std::vector<double>(solution_x_.begin(), solution_x_.end()));
        piecewisebernsteinpoly_result_ = PiecewiseBernsteinPoly(solution_x_2d, tknots, final_time_);
        
        // After calculating final_time_ and bernsteinpoly_result_
        // Flatten bernsteinpoly_result_
        std::vector<double> flattened_result;
        for (const auto& row : piecewisebernsteinpoly_result_) {
            flattened_result.insert(flattened_result.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result, "x.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x_, "x_controlpoints.csv");

        // x2 - g1
        solution_x2_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x2_[i] = g[i];
            //std::cout << "solution_x2_[" << i << "] = " << solution_x2_[i] << std::endl;
        }
        std::vector<std::vector<double>> solution_x2_2d(1, std::vector<double>(solution_x2_.begin(), solution_x2_.end()));      
        piecewisebernsteinpoly_resultx2_ = PiecewiseBernsteinPoly(solution_x2_2d, tknots, final_time_);
        
        std::vector<double> flattened_result1;
        for (const auto& row : piecewisebernsteinpoly_resultx2_) {
            flattened_result1.insert(flattened_result1.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result1, "x1.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x2_, "x1_controlpoints.csv");

        // u - g2
        solution_u_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_u_[i] = g[M_*(N_ + 1) + i];
            //std::cout << "solution_u_[" << i << "] = " << solution_u_[i] << std::endl;
        }
        std::vector<std::vector<double>> solution_u_2d(1, std::vector<double>(solution_u_.begin(), solution_u_.end()));        
        piecewisebernsteinpoly_resultu_ = PiecewiseBernsteinPoly(solution_u_2d, tknots, final_time_);
        
        //std::cout << "solution_u_2d :" << std::endl;
        //for (const auto& row : piecewisebernsteinpoly_resultu_) {
        //    for (const auto& elem : row) {
        //        std::cout << elem << " ";
        //    }
        //    std::cout << std::endl; 
        //}

        std::vector<double> flattened_result2;
        for (const auto& row : piecewisebernsteinpoly_resultu_) {
            flattened_result2.insert(flattened_result2.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result2, "u.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_u_, "u_controlpoints.csv");

        // Save continuity points
        std::vector<double> continuity_times;
        std::vector<double> continuity_values_x;
        std::vector<double> continuity_values_x1;
        std::vector<double> continuity_values_u;

        for (int i = 1; i < M_; ++i) {
            // The time for each continuity point is at the end of each segment
            double time = tknots[i]; // Assuming tknots contain the segment boundaries
            int index = i * (N_ + 1) - 1;

            continuity_times.push_back(time);
            continuity_values_x.push_back(solution_x_[index]);
            continuity_values_x1.push_back(solution_x2_[index]);
            continuity_values_u.push_back(solution_u_[index]);
        }

        writeToCSV(continuity_times, continuity_values_x, "x_continuity.csv");
        writeToCSV(continuity_times, continuity_values_x1, "x1_continuity.csv");
        writeToCSV(continuity_times, continuity_values_u, "u_continuity.csv");

    }
    // Getter for the solution
    const std::vector<Number>& get_solution_x() const { return solution_x_; }

    // Getter for the final objective value
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    int M_;
    double tf_;
    double amax_;
    double amin_;
    double x10_;
    double x1f_;
    double x20_;
    double x2f_;
    PiecewiseBeBOT piecewiseBebot_; 
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
    int N = 50;
    int M = 2;
    double tf = 10;
    double amax = 5.0;
    double amin = -5.0;
    double x10 = -3.0;
    double x1f = 0.0;
    double x20 = 0.0;
    double x2f = 0.0;

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, M, tf, amax, amin, x10, x1f, x20, x2f);
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

// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/pardiso/example_2$ g++ -o pwbebot_example2_v2 ~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/pardiso/example_2/pwbebot_example2.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -I/opt/intel/oneapi/mkl/latest/include -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/pardiso/example_2$ export LD_LIBRARY_PATH=/usr/local/lib/pardiso/panua-pardiso-20230908-linux/lib:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/pardiso/example_2$ ./pwbebot_example2_v2 
