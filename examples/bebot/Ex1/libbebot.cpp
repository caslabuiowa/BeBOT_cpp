#include "../../../../ipopt/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../ipopt/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../include/bebot.h"
#include "../../../include/bernsteinpoly.h"
#include <iomanip>
#include "/opt/intel/oneapi/mkl/latest/include/mkl.h"

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, double amax, double amin, double x10, double x1f, double x20, double x2f)
        : N_(N), tf_(tf), amax_(amax), amin_(amin), x10_(x10), x1f_(x1f), x20_(x20), x2f_(x2f), bebot_(N, tf_) {
        bebot_.calculate();
    }

    void writeToCSV(const std::vector<double>& times, const std::vector<double>& values, const std::string& filename) {
        std::ofstream outFile(filename);

        if (!outFile.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        outFile << "Time,Value\n";
        for (size_t i = 0; i < times.size(); ++i) {
            outFile << std::fixed << std::setprecision(6) << times[i] << "," << values[i] << "\n";
        }

        outFile.close();
    }

    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {
        // Method to allocate the size of the problem 
        // n: number of optimization variables x
        // m: number of constraints g(x)
        // nnz_jac_g: size of nonzero entries of the Jacobian
        // nnz_h_lag: size of nonzero entries of the Hessian
        // index_style: 0-based (C_STYLE) or 1-based (FORTRAN_STYLE)
        n = N_ + 2; // p + t_f
        m = 2 * (N_ + 1); // v + u
        nnz_jac_g = n * m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Method to define the boundary of the variables and constraints
        for (int i = 0; i < N_ + 2; i++) {
            x_l[i] = -std::numeric_limits<double>::infinity();
            x_u[i] = std::numeric_limits<double>::infinity();
        }

        // p(0) = p0
        x_l[0] = x10_;
        x_u[0] = x10_;
        // p(f) = pf
        x_l[(N_ + 2) - 2] = x1f_;
        x_u[(N_ + 2) - 2] = x1f_;
        // 0 <= t <= infinity
        x_l[(N_ + 2) - 1] = 0;
        x_u[(N_ + 2) - 1] = std::numeric_limits<double>::infinity();


        for (int i = 0; i < 2 * (N_ + 1); i++) {
            if (i == 0) { // x2(0) = v0
                g_l[i] = x20_;
                g_u[i] = x20_;
            } else if (i == N_) { // x2(tf) = vf
                g_l[i] = x2f_;
                g_u[i] = x2f_;
            } else if (i >= 1 && i < N_ + 1) {
                g_l[i] = -std::numeric_limits<double>::infinity();
                g_u[i] = std::numeric_limits<double>::infinity();
            } else { // a_min <= u <= a_max
                g_l[i] = amin_;
                g_u[i] = amax_;
            }
        }

        return true;
    }

    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        // Build the initial guess
        for (int i = 0; i < N_ + 2; i++) {
            x[i] = 1; // p(t)
            x[(N_ + 2) - 1] = 10; //tf
        }
        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        // Evaluate objective function
        obj_value = x[(N_ + 2) - 1];
        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Evaluate constraints
        // x2 = x1_dot = x1 * Dm
        Bebot Bebot(N_, x[(N_ + 2) - 1]);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();

        int N_mkl = N_ + 1;
        double* result_flat = new double[N_mkl];

        cblas_dgemv(CblasColMajor, CblasTrans, N_mkl, N_mkl, 1.0, Dm.data(), N_mkl, x, 1, 0.0, result_flat, 1);

        for (int i = 0; i < N_mkl; ++i) {
            g[i] = result_flat[i];
        }

        // u = x2_dot = x2* Dm
        double* result_g2_flat = new double[N_mkl];
        cblas_dgemv(CblasColMajor, CblasTrans, N_mkl, N_mkl, 1.0, Dm.data(), N_mkl, result_flat, 1, 0.0, result_g2_flat, 1);

        for (Index i = 0; i <= N_; i++) {
            g[N_ + 1 + i] = result_g2_flat[i];
        }

        delete[] result_flat;
        delete[] result_g2_flat;

        return true;
    }

    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values) {
        // Evaluate Jacobian of constraints
        if (values == NULL) {
            for (Index i = 0; i < m; i++) {
                for (Index j = 0; j < n; j++) {
                    iRow[i * n + j] = i;
                    jCol[i * n + j] = j;
                }
            }
        }
        return true;
    }

    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
        // Evaluate Jacobian of cost
        return true;
    }

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
        IpoptCalculatedQuantities* ip_cq) { 
        // Extract solution
        // x = optimal solution
        // obj_fun = optimal cost


        // Control points
        solution_x_.resize(n - 1);
        for (Index i = 0; i < n; ++i) {
            solution_x_[i] = x[i];
        }
        final_obj_value_ = obj_value;

        // Print the solution
        std::cout << "Optimal Solution: ";
        for (int i = 0; i < n - 1; ++i) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;

        // Retrieve the final objective value
        std::cout << "Final Objective Value: " << obj_value << std::endl;



        // Tf
        tf_ = x[(N_ + 2) - 1];
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();
        
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * obj_value / 999.0;
        }
        
        std::vector<std::vector<double>> solution_x_2d(1, std::vector<double>(solution_x_.begin(), solution_x_.end()));
        bernsteinpoly_result_ = BernsteinPoly(solution_x_2d, final_time_, 0, tf_);
        
        std::vector<double> flattened_result;
        for (const auto& row : bernsteinpoly_result_) {
            flattened_result.insert(flattened_result.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result, "x.csv");
        writeToCSV(bebot_.getNodes(), solution_x_, "x_controlpoints.csv");

        solution_x2_.resize(n - 1);
        for (Index i = 0; i < n - 1; ++i) {
            solution_x2_[i] = g[i];
        }
        std::vector<std::vector<double>> solution_x2_2d(1, std::vector<double>(solution_x2_.begin(), solution_x2_.end()));      
        bernsteinpoly_resultx2_ = BernsteinPoly(solution_x2_2d, final_time_, 0, tf_);
        
        std::vector<double> flattened_result1;
        for (const auto& row : bernsteinpoly_resultx2_) {
            flattened_result1.insert(flattened_result1.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result1, "x1.csv");
        writeToCSV(bebot_.getNodes(), solution_x2_, "x1_controlpoints.csv");

        solution_u_.resize(n - 1);
        for (Index i = 0; i < n - 1; ++i) {
            solution_u_[i] = g[N_ + 1 + i];
        }
        std::vector<std::vector<double>> solution_u_2d(1, std::vector<double>(solution_u_.begin(), solution_u_.end()));        
        bernsteinpoly_resultu_ = BernsteinPoly(solution_u_2d, final_time_, 0, tf_);
        
        std::vector<double> flattened_result2;
        for (const auto& row : bernsteinpoly_resultu_) {
            flattened_result2.insert(flattened_result2.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result2, "u.csv");
        writeToCSV(bebot_.getNodes(), solution_u_, "u_controlpoints.csv");
    }

    const std::vector<Number>& get_solution_x() const { return solution_x_; }
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

extern "C" {
    PointSetProblem* create_point_set_problem(int N, double tf, double amax, double amin, double x10, double x1f, double x20, double x2f) {
        std::cout << "Received parameters:" << std::endl;
        std::cout << "N: " << N << std::endl;
        std::cout << "tf: " << tf << std::endl;
        std::cout << "amax: " << amax << std::endl;
        std::cout << "amin: " << amin << std::endl;
        std::cout << "x10: " << x10 << std::endl;
        std::cout << "x1f: " << x1f << std::endl;
        std::cout << "x20: " << x20 << std::endl;
        std::cout << "x2f: " << x2f << std::endl;
        return new PointSetProblem(N, tf, amax, amin, x10, x1f, x20, x2f);
    }

    void solve_point_set_problem(PointSetProblem* problem) {
        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetStringValue("linear_solver", "ma57");
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
        app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        app->Options()->SetIntegerValue("max_iter", 5000);
        app->Options()->SetNumericValue("tol", 1e-6);

        app->RethrowNonIpoptException(true);
        ApplicationReturnStatus status = app->Initialize();
        if (status != Solve_Succeeded) {
            std::cerr << "IPOPT initialization failed!" << std::endl;
            return;
        }

        status = app->OptimizeTNLP(problem);

        if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
            std::cout << "Optimization succeeded!" << std::endl;
        } else {
            std::cerr << "Optimization failed with status " << status << std::endl;
        }
    }

    void get_solution(PointSetProblem* problem, double* solution, int n) {
        const std::vector<double>& sol = problem->get_solution_x();
        std::copy(sol.begin(), sol.end(), solution);
    }

    double get_final_objective_value(PointSetProblem* problem) {
        return problem->get_final_obj_value();
    }

    void destroy_point_set_problem(PointSetProblem* problem) {
        delete problem;
    }
}
