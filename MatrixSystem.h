#ifndef MATRIXSYSTEM_H_INCLUDED
#define MATRIXSYSTEM_H_INCLUDED

#include "eigen-3.4.0/Eigen/Dense"
#include <vector>
#include "BoundaryCondition.h"
#include "CoefEquation.h"

namespace system_matrix
{
    class MatrixSystem
    {
    public:
        MatrixSystem(coef_eq::CoefEquation alpha,
                     coef_eq::CoefEquation beta,
                     coef_eq::CoefEquation gamma,
                     coef_eq::CoefEquation delta,
                     double theta,
                     boundary::BoundaryCondition boundary_small_spot,
                     boundary::BoundaryCondition boundary_big_spot,
                     std::vector<double> Xt1);
        ~MatrixSystem() = default;
        void solve();
        std::vector<double> get_result();

    private:
        Eigen::MatrixXd m_A_prime;
        Eigen::MatrixXd m_A_second;
        Eigen::MatrixXd m_b;
        std::vector<double> m_result; // use of vector in order to avoid having to import eigen in client

        double Omega(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double a_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double b_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double c_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double d_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double e_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;

        double Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;

        double Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double a_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;

    };
}


#endif // MATRIXSYSTEM_H_INCLUDED
