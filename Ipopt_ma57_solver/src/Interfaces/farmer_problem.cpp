#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"


using namespace Ipopt;

class FarmerProblem : public TNLP {
public:
    // Constructor
    FarmerProblem(double fence_length) : fence_length_(fence_length) {}

    // Define problem dimensions
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {
        // Number of variables (x1, x2)
        n = 2;
        
        // Number of constraints (1 for the fence length)
        m = 1;
        
        // Number of non-zero elements in the Jacobian of the constraints
        nnz_jac_g = 2;
        
        // Number of non-zero elements in the Hessian of the Lagrangian (0 if not available)
        nnz_h_lag = 0;
        
        // Use C-style indexing
        index_style = TNLP::C_STYLE;
        return true;
    }

    // Define bounds for variables and constraints
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Variable bounds
        x_l[0] = 0.0;  // x1 >= 0
        x_l[1] = 0.0;  // x2 >= 0
        x_u[0] = fence_length_;  // x1 <= a
        x_u[1] = fence_length_ / 2.0;  // x2 <= a/2

        // Constraint bounds
        g_l[0] = 0.0;  // g(x1, x2) >= 0
        g_u[0] = 0.0;  // g(x1, x2) <= 0

        return true;
    }

    // Define initial guess for the variables
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        x[0] = 0.0;  // Initial guess for x1
        x[1] = 0.0;  // Initial guess for x2
        return true;
    }

    // Define the objective function
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        obj_value = -x[0] * x[1];  // J(x, a) = -x1 * x2
        return true;
    }

    // Define the constraints
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        g[0] = x[0] + 2.0 * x[1] - fence_length_;  // x1 + 2*x2 - a
        return true;
    }

    // Define the Jacobian of the constraints
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values) {
        if (values == NULL) {
            // Structure of the Jacobian
            iRow[0] = 0;  // Row index for the constraint
            jCol[0] = 0;  // Column index for x1
            iRow[1] = 0;  // Row index for the constraint
            jCol[1] = 1;  // Column index for x2
        } else {
            // Values of the Jacobian
            values[0] = 1.0;  // Partial derivative of g wrt x1
            values[1] = 2.0;  // Partial derivative of g wrt x2
        }
        return true;
    }
private:
    double fence_length_;
};

int main() {
    // Fixed fence length
    double a = 1000.0;

    // Create an instance of the IPOPT application
    SmartPtr<TNLP> farmerProblem = new FarmerProblem(a);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // Initialize IPOPT and set options
    app->RethrowNonIpoptException(true);
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << "IPOPT initialization failed!" << std::endl;
        return -1;
    }

    // Solve the optimization problem
    status = app->OptimizeTNLP(farmerProblem);

    if (status == Solve_Succeeded) {
        // Retrieve and print the solution
        Number obj_value;
        const Number* x = farmerProblem->x_sol();
        farmerProblem->eval_f(2, x, true, obj_value);

        std::cout << "Optimal x1: " << x[0] << std::endl;
        std::cout << "Optimal x2: " << x[1] << std::endl;
        std::cout << "Optimal J(x, a): " << obj_value << std::endl;
    } else {
        std::cout << "IPOPT optimization failed with status " << status << std::endl;
    }

    return 0;
}
