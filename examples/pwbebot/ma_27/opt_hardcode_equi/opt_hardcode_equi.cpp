#include "../../../../Ipopt_ma27_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma27_solver/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/piecewisebebot.h"
#include "../../../../include/piecewisebernsteinpoly.h"
#include "../../../../include/piecewiseintegrationmatrix.h"
#include "../../../../include/piecewisebernsteinproduct.h"
#include "../../../../include/knots.h"
#include <iomanip>
#include "mkl.h"


using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int K, double t0, double tf, double amax, double amin, double x10, double x1f, double x20, double x2f)
        : K_(K), t0_(t0), tf_(tf), amax_(amax), amin_(amin), x10_(x10), x1f_(x1f), x20_(x20), x2f_(x2f), piecewiseBebot_(1, initializeTknots()) {
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
        n = (3 * K_) + 5; 
        m = (12 * K_) + 5 + 2;
        nnz_jac_g = n*m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Define variable bounds
        for (Index i = 0; i < (3 * K_) + 5; i++) {
            /*if (i == 0) {
                x_l[i] = x10_;
                x_u[i] = x10_;
            } else if (i == K_ - 1) {
                x_l[i] = x1f_;
                x_u[i] = x1f_;
            } else if (i == K_) {
                x_l[i] = x20_;
                x_u[i] = x20_;
            } else if (i == 2 * K_ - 1) {
                x_l[i] = x2f_;
                x_u[i] = x2f_;
            } else if (i >= 2 * K_ && i < 3 * K_) {
                x_l[i] = amin_;
                x_u[i] = amax_;
            } else {*/
                x_l[i] = -std::numeric_limits<double>::infinity();
                x_u[i] = std::numeric_limits<double>::infinity();
            //}
            //std::cout << "x_l[" << i << "] = " << x_u[i] << std::endl;
            //std::cout << "x_u[" << i << "] = " << x_u[i] << std::endl;
            
        }

        // Print the bounds for each variable
        //for (Index i = 0; i < M_ * (N_ + 1) + 1; ++i) {
        //    std::cout << "x_l[" << i << "] = " << x_l[i] << ", x_u[" << i << "] = " << x_u[i] << std::endl;
        //}

        // Define constraint bounds
        for (Index i = 0; i < (12 * K_) + 5 + 2; i++) {
            
            if (i == 0) {
                g_l[i] = x10_;
                g_u[i] = x10_;
            } else if (i == 3 * K_ - 1) {
                g_l[i] = x1f_;
                g_u[i] = x1f_;
            } else if (i == 3 * K_) {
                g_l[i] = x20_;
                g_u[i] = x20_;
            } else if (i == 6 * K_ - 1) {
                g_l[i] = x2f_;
                g_u[i] = x2f_;
            } else if (i >= 6 * K_ && i < 8 * K_ - 1) {
                g_l[i] = amin_;
                g_u[i] = amax_;
            } else if (i == 8 * K_ && i < 10 * K_ - 1) {
                g_l[i] = 0;
                g_u[i] = 0;
            } else if (i == (12 * K_) + 5 || i == (12 * K_) + 5 + 1) {
                g_l[i] = 0;
                g_u[i] = 0;
            } else {
                g_l[i] = -std::numeric_limits<double>::infinity();
                g_u[i] = std::numeric_limits<double>::infinity();
            }
            
            std::cout << "g_l[" << i << "] = " << g_u[i] << std::endl;
            std::cout << "g_u[" << i << "] = " << g_u[i] << std::endl;
        }


        // Print constraint bounds
        //for (Index i = 0; i < M_ * 2 * (N_ + 1) + M_-1; ++i) {
        //    std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        //}

        return true;
    }

    // initialization of the starting point
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (int i = 0; i < (3 * K_) + 5; i++) {
            x[i] = 1;
        }

        //for (Index i = 0; i < M_ * (N_ + 1) + 1; ++i) {
        //    std::cout << "x[" << i << "] = " << x[i] << std::endl;
        //}

        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        // Initialize tknots_ directly here
        tknots_.resize(K_ + 1);  // Resize to hold K+1 elements
        for (int i = 0; i <= K_; ++i) {  // Loop from 0 to K inclusive
            tknots_[i] = t0_ + i * (tf_ - t0_) / K_;  // Divide by K to include the endpoint
            //std::cout << "tknots_[" << i << "] = " << tknots_[i] << std::endl;
        }
        // Create integration matrices I_1 and I_2
        std::vector<std::vector<double>> I_1 = PiecewiseIntegrationMatrix(K_, 0, tknots_);
        std::vector<std::vector<double>> I_2 = PiecewiseIntegrationMatrix(K_, 1, tknots_);

        // Print I_1
        //std::cout << "I_1 Matrix in eval_g:" << std::endl;
        //for (const auto& row : I_1) {
        //    for (const auto& elem : row) {
        //        std::cout << std::setw(10) << std::setprecision(4) << elem << " ";
        //    }
        //    std::cout << std::endl;
        //}

        // Print I_2
        //std::cout << "I_2 Matrix in eval_g:" << std::endl;
        //for (const auto& row : I_2) {
        //    for (const auto& elem : row) {
        //        std::cout << std::setw(10) << std::setprecision(4) << elem << " ";
        //    }
        //    std::cout << std::endl;
        //}

        // Extract vectors from x
        std::vector<double> x1_z(x, x + K_);  // vector x from 0 to K-1
        std::vector<double> x2_z(x + K_, x + 2 * K_);  // vector x from K to 2*K-1
        std::vector<double> u_z(x + 2 * K_, x + 3 * K_);  // vector x from 2*K to 3*K-1
        /*
        // Print x1_z vector
        std::cout << "x1_z values: ";
        for (auto val : x1_z) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print x2_z vector
        std::cout << "x2_z values: ";
        for (auto val : x2_z) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print u_z vector
        std::cout << "u_z values: ";
        for (auto val : u_z) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        double x1pc = x[3 * K_];
        double x1c = x[3 * K_ + 1];
        double x2pc = x[3 * K_ + 2];
        double x2c = x[3 * K_ + 3];
        double uc = x[3 * K_ + 4];

        //std::cout << "x1pc: " << x1pc << std::endl;
        //std::cout << "x1c: " << x1c << std::endl;
        //std::cout << "x2pc: " << x2pc << std::endl;
        //std::cout << "x2c: " << x2c << std::endl;
        //std::cout << "uc: " << uc << std::endl;

        int rows_I1 = I_1.size();        // Number of rows in I_1
        int cols_I1 = I_1[0].size();     // Number of columns in I_1
        int rows_I2 = I_2.size();        // Number of rows in I_2
        int cols_I2 = I_2[0].size();     // Number of columns in I_2

        // Initialize vectors based on matrix dimensions
        std::vector<double> dx1(cols_I1, 0.0), dx2(cols_I1, 0.0), u(cols_I1, 0.0);
        std::vector<double> x1(cols_I2, 0.0), x2(cols_I2, 0.0);
        
        std::vector<double> flat_I1(I_1.size() * I_1[0].size());
        for (int i = 0; i < I_1.size(); ++i) {
            std::copy(I_1[i].begin(), I_1[i].end(), flat_I1.begin() + i * I_1[0].size());
        }

        std::vector<double> flat_I2(I_2.size() * I_2[0].size());
        for (int i = 0; i < I_2.size(); ++i) {
            std::copy(I_2[i].begin(), I_2[i].end(), flat_I2.begin() + i * I_2[0].size());
        }

        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, flat_I1.data(), I_1[0].size(), x1_z.data(), 1, 0, dx1.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, flat_I1.data(), I_1[0].size(), x2_z.data(), 1, 0, dx2.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, flat_I1.data(), I_1[0].size(), u_z.data(), 1, 0, u.data(), 1);

        /*
        // Print dx1 vector after computation
        std::cout << "dx1 values after computation: ";
        for (auto val : dx1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print dx2 vector after computation
        std::cout << "dx2 values after computation: ";
        for (auto val : dx2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print u vector after computation
        std::cout << "u values after computation: ";
        for (auto val : u) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */

        // Create a vector of ones for use with daxpy
        std::vector<double> ones(dx1.size(), 1.0);

        cblas_daxpy(dx1.size(), x1pc, ones.data(), 1, dx1.data(), 1);
        cblas_daxpy(dx2.size(), x2pc, ones.data(), 1, dx2.data(), 1);
        cblas_daxpy(u.size(), uc, ones.data(), 1, u.data(), 1);

        /*

        // Print dx1 vector after addition
        std::cout << "dx1 values after addition: ";
        for (auto val : dx1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print dx2 vector after addition
        std::cout << "dx2 values after addition: ";
        for (auto val : dx2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print u vector after addition
        std::cout << "u values after addition: ";
        for (auto val : u) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        */

        cblas_dgemv(CblasRowMajor, CblasTrans, I_2.size(), I_2[0].size(), 1, flat_I2.data(), I_2[0].size(), dx1.data(), 1, 0, x1.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_2.size(), I_2[0].size(), 1, flat_I2.data(), I_2[0].size(), dx2.data(), 1, 0, x2.data(), 1);

        /*

        // Print x1 vector after computation
        std::cout << "x1 values after computation: ";
        for (auto val : x1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print x2 vector after computation
        std::cout << "x2 values after computation: ";
        for (auto val : x2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */

        std::vector<double> singles(x1.size(), 1.0);

        cblas_daxpy(x1.size(), x1c, singles.data(), 1, x1.data(), 1);
        cblas_daxpy(x2.size(), x2c, singles.data(), 1, x2.data(), 1);

        /*
        // Print dx1 vector after addition
        std::cout << "x1 values after addition: ";
        for (auto val : x1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print dx2 vector after addition
        std::cout << "x2 values after addition: ";
        for (auto val : x2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        //for (auto& val : x1) val += x1c;
        //for (auto& val : x2) val += x2c;

        std::vector<double> u_product = PiecewiseBernsteinProduct(u, u, K_, 1);
        
        /*
        std::cout << "u_product values: ";
        for (auto val : u_product) {
            std::cout << std::setprecision(4) << val << " ";
        }
        std::cout << std::endl;
        */
        
        obj_value = cblas_dasum(u_product.size(), &u_product[0], 1);
        
        //std::cout << "Sum of absolute values in u_product: " << std::setprecision(4) << obj_value << std::endl;
        //std::cout << "objective value = " << obj_value << std::endl;
        return true;
    }
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    
    //G vector/*
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Initialize tknots_ directly here
        tknots_.resize(K_ + 1);  // Resize to hold K+1 elements
        for (int i = 0; i <= K_; ++i) {  // Loop from 0 to K inclusive
            tknots_[i] = t0_ + i * (tf_ - t0_) / K_;  // Divide by K to include the endpoint
            //std::cout << "tknots_[" << i << "] = " << tknots_[i] << std::endl;
        }
        // Create integration matrices I_1 and I_2
        std::vector<std::vector<double>> I_1 = PiecewiseIntegrationMatrix(K_, 0, tknots_);
        std::vector<std::vector<double>> I_2 = PiecewiseIntegrationMatrix(K_, 1, tknots_);

        // Extract vectors from x
        std::vector<double> gx1_z(x, x + K_);  // vector x from 0 to K-1
        std::vector<double> gx2_z(x + K_, x + 2 * K_);  // vector x from K to 2*K-1
        std::vector<double> gu_z(x + 2 * K_, x + 3 * K_);  // vector x from 2*K to 3*K-1

        double gx1pc = x[3 * K_];
        double gx1c = x[3 * K_ + 1];
        double gx2pc = x[3 * K_ + 2];
        double gx2c = x[3 * K_ + 3];
        double guc = x[3 * K_ + 4];


        int rows_I1 = I_1.size();        // Number of rows in I_1
        int cols_I1 = I_1[0].size();     // Number of columns in I_1
        int rows_I2 = I_2.size();        // Number of rows in I_2
        int cols_I2 = I_2[0].size();     // Number of columns in I_2

        // Initialize vectors based on matrix dimensions
        std::vector<double> gdx1(cols_I1, 0.0), gdx2(cols_I1, 0.0), gu(cols_I1, 0.0); 
        std::vector<double> gx1(cols_I2, 0.0), gx2(cols_I2, 0.0);
        
        std::vector<double> gflat_I1(I_1.size() * I_1[0].size());
        for (int i = 0; i < I_1.size(); ++i) {
            std::copy(I_1[i].begin(), I_1[i].end(), gflat_I1.begin() + i * I_1[0].size());
        }

        std::vector<double> gflat_I2(I_2.size() * I_2[0].size());
        for (int i = 0; i < I_2.size(); ++i) {
            std::copy(I_2[i].begin(), I_2[i].end(), gflat_I2.begin() + i * I_2[0].size());
        }

        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, gflat_I1.data(), I_1[0].size(), gx1_z.data(), 1, 0, gdx1.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, gflat_I1.data(), I_1[0].size(), gx2_z.data(), 1, 0, gdx2.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, gflat_I1.data(), I_1[0].size(), gu_z.data(), 1, 0, gu.data(), 1);

        // Create a vector of ones for use with daxpy
        std::vector<double> gones(gdx1.size(), 1.0);

        cblas_daxpy(gdx1.size(), gx1pc, gones.data(), 1, gdx1.data(), 1);
        cblas_daxpy(gdx2.size(), gx2pc, gones.data(), 1, gdx2.data(), 1);
        cblas_daxpy(gu.size(), guc, gones.data(), 1, gu.data(), 1);
        /*
        std::cout << "gdx1: ";
        for (auto val : gdx1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "gdx2: ";
        for (auto val : gdx2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print u vector after addition
        std::cout << "gu: ";
        for (auto val : gu) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */

        cblas_dgemv(CblasRowMajor, CblasTrans, I_2.size(), I_2[0].size(), 1, gflat_I2.data(), I_2[0].size(), gdx1.data(), 1, 0, gx1.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_2.size(), I_2[0].size(), 1, gflat_I2.data(), I_2[0].size(), gdx2.data(), 1, 0, gx2.data(), 1);

        std::vector<double> gsingles(gx1.size(), 1.0);

        cblas_daxpy(gx1.size(), gx1c, gsingles.data(), 1, gx1.data(), 1);
        cblas_daxpy(gx2.size(), gx2c, gsingles.data(), 1, gx2.data(), 1);
        /*
        std::cout << "gx1: ";
        for (auto val : gx1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "gx2: ";
        for (auto val : gx2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        PiecewiseBeBOT piecewiseBebot(1, tknots_);
        piecewiseBebot.calculate();

        std::vector<std::vector<double>> Dm_2D = piecewiseBebot.getDifferentiationMatrix();
        
        // Print the 2D Differentiation Matrix
        //std::cout << "2D Differentiation Matrix (Dm_2D):" << std::endl;
        //for (const auto& row : Dm_2D) {
        //    for (const auto& val : row) {
        //        std::cout << std::fixed << std::setprecision(4) << val << " ";
        //    }
        //    std::cout << std::endl;
        //}

        // Flatten the 2D matrix into a 1D vector
        std::vector<double> Dm_flat;
        int numRows = Dm_2D.size();
        int numCols = numRows > 0 ? Dm_2D[0].size() : 0;
        Dm_flat.reserve(numRows * numCols);

        for (const auto& row : Dm_2D) {
            Dm_flat.insert(Dm_flat.end(), row.begin(), row.end());
        }

        //std::cout << "Flattened Differentiation Matrix (Dm_flat):" << std::endl;
        //for (double val : Dm_flat) {
        //    std::cout << std::fixed << std::setprecision(4) << val << " ";
        //}
        //std::cout << std::endl;


        std::vector<double> gu_dot(numRows - 1, 0.0);

        // Perform the matrix-vector multiplication
        cblas_dgemv(CblasRowMajor, CblasTrans, numRows, numCols, 1.0, Dm_flat.data(), numCols, gu.data(), 1, 0.0, gu_dot.data(), 1);

        // Print u vector after addition
        //std::cout << "gu_dot: ";
        //for (auto val : gu_dot) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        // Getting Knots
        std::vector<double> knots_gdx1 = Knots(gdx1, K_);
        std::vector<double> knots_gx2 = Knots(gx2, K_);
        std::vector<double> knots_gdx2 = Knots(gdx2, K_);
        std::vector<double> knots_gu = Knots(gu, K_);

        // Print knots
        //std::cout << "knots_gdx1: ";
        //for (auto val : knots_gdx1) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        //std::cout << "knots_gx2: ";
        //for (auto val : knots_gx2) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        //std::cout << "knots_gdx2: ";
        //for (auto val : knots_gdx2) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        //std::cout << "knots_gu: ";
        //for (auto val : knots_gu) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;


        std::vector<Number> g_eq1(2 * K_, 0.0);
        std::vector<Number> g_eq2(2 * K_, 0.0);

        // Calculate g_eq1 = knots_gdx1 - knots_gx2
        vdSub(2 * K_, knots_gdx1.data(), knots_gx2.data(), g_eq1.data());

        // Calculate g_eq2 = knots_gdx2 - knots_gu
        vdSub(2 * K_, knots_gdx2.data(), knots_gu.data(), g_eq2.data());
        
        //std::cout << "g_eq1: ";
        //for (auto val : g_eq1) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;
    
        //std::cout << "g_eq2: ";
        //for (auto val : g_eq2) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        

        //std::cout << "Complete g vector:" << std::endl;
        //for (Index i = 0; i < M_ * 2 * (N_ + 1) +M_-1; ++i) {
        //    std::cout << "g[" << i << "] = " << g[i] << std::endl;
        //}

        // x1 vectors
        for (Index i = 0; i < 3 * K_; ++i) {
            g[i] = gx1[i];  // Assuming 'x1' is defined and populated elsewhere
            std::cout << "g[" << i << "] = " << g[i] << std::endl;
            
        }

        // x2 vectors
        for (Index i = 0; i < 3 * K_; ++i) {
            g[3 * K_ + i] = gx2[i];  // Assuming 'x2' is defined and populated elsewhere
            std::cout << "g[" << 3 * K_ + i << "] = " << g[3 * K_ + i] << std::endl;
            
        }
        // u vector
        for (Index i = 0; i < 2 * K_; ++i) {
            g[6 * K_ + i] = gu[i];  
            std::cout << "g[" << 6 * K_ + i << "] = " << g[6 * K_ + i] << std::endl;     
        }

        // knots difference vectors
        for (Index i = 0; i < 2 * K_; ++i) {
            g[8 * K_ + i] = g_eq1[i];
            std::cout << "g[" << 8 * K_ + i << "] = " << g[8 * K_ + i] << std::endl;  
        }

        // knots difference vectors
        for (Index i = 0; i < 2 * K_; ++i) {
            g[10 * K_ + i] = g_eq2[i];
            std::cout << "g[" << 10 * K_ + i << "] = " << g[10 * K_ + i] << std::endl;  
        }
        
        // constants
        g[(12 * K_)] = gx1pc;
        g[(12 * K_) + 1] = gx1c;
        g[(12 * K_) + 2] = gx2pc;
        g[(12 * K_) + 3] = gx2c;
        g[(12 * K_) + 4] = guc;

        // u_dot difference udot(2) - udot(3);
        g[(12 * K_) + 5] = g[(6 * K_) + 1] - g[(6 * K_) + 2];

        // u_dot difference udot(end-2) - udot(end-1);
        g[(12 * K_) + 5 + 1] = g[(8 * K_) - 3] - g[(8 * K_) - 2];

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
    {   final_obj_value_ = obj_value;
        
        // Initialize tknots_ directly here
        tknots_.resize(K_ + 1);  // Resize to hold K+1 elements
        for (int i = 0; i <= K_; ++i) {  // Loop from 0 to K inclusive
            tknots_[i] = t0_ + i * (tf_ - t0_) / K_;  // Divide by K to include the endpoint
            //std::cout << "tknots_[" << i << "] = " << tknots_[i] << std::endl;
        }
        // Create integration matrices I_1 and I_2
        std::vector<std::vector<double>> I_1 = PiecewiseIntegrationMatrix(K_, 0, tknots_);
        std::vector<std::vector<double>> I_2 = PiecewiseIntegrationMatrix(K_, 1, tknots_);


        // Extract vectors from x
        std::vector<double> x1_z(x, x + K_);  // vector x from 0 to K-1
        std::vector<double> x2_z(x + K_, x + 2 * K_);  // vector x from K to 2*K-1
        std::vector<double> u_z(x + 2 * K_, x + 3 * K_);  // vector x from 2*K to 3*K-1

        double x1pc = x[3 * K_];
        double x1c = x[3 * K_ + 1];
        double x2pc = x[3 * K_ + 2];
        double x2c = x[3 * K_ + 3];
        double uc = x[3 * K_ + 4];

        int rows_I1 = I_1.size();        // Number of rows in I_1
        int cols_I1 = I_1[0].size();     // Number of columns in I_1
        int rows_I2 = I_2.size();        // Number of rows in I_2
        int cols_I2 = I_2[0].size();     // Number of columns in I_2

        // Initialize vectors based on matrix dimensions
        std::vector<double> dx1(cols_I1, 0.0), dx2(cols_I1, 0.0), u(cols_I1, 0.0);
        std::vector<double> x1(cols_I2, 0.0), x2(cols_I2, 0.0);
        
        std::vector<double> flat_I1(I_1.size() * I_1[0].size());
        for (int i = 0; i < I_1.size(); ++i) {
            std::copy(I_1[i].begin(), I_1[i].end(), flat_I1.begin() + i * I_1[0].size());
        }

        std::vector<double> flat_I2(I_2.size() * I_2[0].size());
        for (int i = 0; i < I_2.size(); ++i) {
            std::copy(I_2[i].begin(), I_2[i].end(), flat_I2.begin() + i * I_2[0].size());
        }

        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, flat_I1.data(), I_1[0].size(), x1_z.data(), 1, 0, dx1.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, flat_I1.data(), I_1[0].size(), x2_z.data(), 1, 0, dx2.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_1.size(), I_1[0].size(), 1, flat_I1.data(), I_1[0].size(), u_z.data(), 1, 0, u.data(), 1);

        std::vector<double> ones(dx1.size(), 1.0);

        cblas_daxpy(dx1.size(), x1pc, ones.data(), 1, dx1.data(), 1);
        cblas_daxpy(dx2.size(), x2pc, ones.data(), 1, dx2.data(), 1);
        cblas_daxpy(u.size(), uc, ones.data(), 1, u.data(), 1);

        cblas_dgemv(CblasRowMajor, CblasTrans, I_2.size(), I_2[0].size(), 1, flat_I2.data(), I_2[0].size(), dx1.data(), 1, 0, x1.data(), 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, I_2.size(), I_2[0].size(), 1, flat_I2.data(), I_2[0].size(), dx2.data(), 1, 0, x2.data(), 1);

        std::vector<double> singles(x1.size(), 1.0);

        cblas_daxpy(x1.size(), x1c, singles.data(), 1, x1.data(), 1);
        cblas_daxpy(x2.size(), x2c, singles.data(), 1, x2.data(), 1);

        PiecewiseBeBOT piecewiseBebot(1, tknots_);
        piecewiseBebot.calculate();

        std::vector<std::vector<double>> Dm_2D = piecewiseBebot.getDifferentiationMatrix();
        
        // Flatten the 2D matrix into a 1D vector
        std::vector<double> Dm_flat;
        int numRows = Dm_2D.size();
        int numCols = numRows > 0 ? Dm_2D[0].size() : 0;
        Dm_flat.reserve(numRows * numCols);

        for (const auto& row : Dm_2D) {
            Dm_flat.insert(Dm_flat.end(), row.begin(), row.end());
        }

        std::vector<double> u_dot(numRows - 1, 0.0);

        // Perform the matrix-vector multiplication
        cblas_dgemv(CblasRowMajor, CblasTrans, numRows, numCols, 1.0, Dm_flat.data(), numCols, u.data(), 1, 0.0, u_dot.data(), 1);

        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * (tf_ - t0_) / 999.0;
            //std::cout << "final_time_[" << i << "] = " << final_time_[i] << std::endl;
        }
        
        // x1
        std::vector<std::vector<double>> x1_(1, std::vector<double>(x1.begin(), x1.end()));
        piecewisebernsteinpoly_result_ = PiecewiseBernsteinPoly(x1_, tknots_, final_time_);
        
        
        // Flatten bernsteinpoly_result_
        std::vector<double> flattened_result;
        for (const auto& row : piecewisebernsteinpoly_result_) {
            flattened_result.insert(flattened_result.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result, "x1.csv");
        writeToCSV(tknots_, x1, "x1_controlpoints.csv");

        
        // x2
        std::vector<std::vector<double>> x2_(1, std::vector<double>(x2.begin(), x2.end()));
        piecewisebernsteinpoly_resultx2_ = PiecewiseBernsteinPoly(x2_, tknots_, final_time_);
                
        // Flatten bernsteinpoly_result_
        std::vector<double> flattened_result1;
        for (const auto& row : piecewisebernsteinpoly_resultx2_) {
            flattened_result1.insert(flattened_result1.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result1, "x2.csv");
        writeToCSV(piecewiseBebot_.getNodes(), x2, "x2_controlpoints.csv");


        // x2
        std::vector<std::vector<double>> u_(1, std::vector<double>(u.begin(), u.end()));
        piecewisebernsteinpoly_resultu_ = PiecewiseBernsteinPoly(u_, tknots_, final_time_);
                
        // Flatten bernsteinpoly_result_
        std::vector<double> flattened_resultu;
        for (const auto& row : piecewisebernsteinpoly_resultu_) {
            flattened_resultu.insert(flattened_resultu.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_resultu, "u.csv");
        writeToCSV(piecewiseBebot_.getNodes(), u, "u_controlpoints.csv");

        


    }
    // Getter for the solution
    const std::vector<Number>& get_solution_x() const { return solution_x_; }

    // Getter for the final objective value
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int K_;
    double t0_;
    double tf_;
    double amax_;
    double amin_;
    double x10_;
    double x1f_;
    double x20_;
    double x2f_;
    std::vector<double> tknots_;

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
    std::vector<double> initializeTknots() {
        std::vector<double> tknots(K_);
        for (int i = 0; i < K_; ++i) {
            tknots[i] = t0_ + i * (tf_ - t0_) / (K_ - 1);
            std::cout << "tknots_[" << i << "] = " << tknots[i] << std::endl;
        }
        return tknots;
    }

public:
    
};

int main() {
    int K = 5;
    double t0 = 0;
    double tf = 10;
    double amax = 5.0;
    double amin = -5.0;
    double x10 = -3.0;
    double x1f = 0.0;
    double x20 = 0.0;
    double x2f = 0.0;

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(K, t0, tf, amax, amin, x10, x1f, x20, x2f);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    //app->Options()->SetStringValue("linear_solver", "ma27");
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

// g++ -o opt_hardcode_equi opt_hardcode_equi.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewiseintegrationmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebernsteinproduct.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/knots.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinproduct.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I../../../../Ipopt_ma27_solver/include/coin-or -I/opt/intel/oneapi/mkl/latest/include -L../../../../Ipopt_ma27_solver/lib -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lipopt -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++