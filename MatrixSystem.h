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
        MatrixSystem(coef_eq::CoefEquation *alpha, coef_eq::CoefEquation *beta, coef_eq::CoefEquation* gamma, coef_eq::CoefEquation* delta,
                     double theta, double dt, double dx, double sigma, double r,
                     boundary::BoundaryCondition* boundary_small_spot, boundary::BoundaryCondition* boundary_big_spot,
                     std::vector<double> Xt1, double spot_min, double spot_max, double S0, double loop, double maturity);
        ~MatrixSystem() = default;
        std::vector<double> solve();

    private:
        //Omega*Xt = A'*Xt+1 + A''*Xt +b <==> m_A*Xt = m_b
        Eigen::MatrixXd m_A; //avec m_A = (Omega - A'')
        Eigen::MatrixXd m_b; //et m_b = A'*Xt+1 + b

        // not used yet but can be useful if coeff equations become variables
        double m_S0;
        double m_loop;
        double m_maturity;
        double m_spot_min;
        double m_spot_max;

        double Omega(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double a_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double b_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double c_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double d_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double e_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
    };
}


#endif // MATRIXSYSTEM_H_INCLUDED
