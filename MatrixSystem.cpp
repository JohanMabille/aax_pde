#include "MatrixSystem.h"
#include "CoefEquation.h"
#include "BoundaryCondition.h"
#include <vector>
#include <cmath>
#include "eigen-3.4.0/Eigen/Dense"

namespace system_matrix
{
    MatrixSystem::MatrixSystem(coef_eq::CoefEquation alpha,
                     coef_eq::CoefEquation beta,
                     coef_eq::CoefEquation gamma,
                     coef_eq::CoefEquation delta,
                     double theta,
                     boundary::BoundaryCondition boundary_small_spot,
                     boundary::BoundaryCondition boundary_big_spot,
                     std::vector<double> Xt1)
    {
        //todo
    }

    std::vector<double> MatrixSystem::get_result()
    {
        return m_result;
    }


    //i != 0, i!=N

    double MatrixSystem::Omega(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = theta * dt * (2*alpha.get_value(args) / pow(dx, 2.0) - gamma.get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::a_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (-2*alpha.get_value(args) / pow(dx, 2.0) + gamma.get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::b_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = theta * dt * (alpha.get_value(args) / pow(dx, 2.0) - beta.get_value(args) / 2*dx);
        return res;
    }

    double MatrixSystem::c_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = theta * dt * (alpha.get_value(args) / pow(dx, 2.0) + beta.get_value(args) / 2*dx);
        return res;
    }

    double MatrixSystem::d_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (alpha.get_value(args) / pow(dx, 2.0) + beta.get_value(args) / 2*dx);
        return res;
    }

    double MatrixSystem::e_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (alpha.get_value(args) / pow(dx, 2.0) - beta.get_value(args) / 2*dx);
        return res;
    }

    //i = 0

    double MatrixSystem::Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = -theta * dt * (alpha.get_value(args) / pow(dx, 2.0) - beta.get_value(args) / dx + gamma.get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (alpha.get_value(args) / pow(dx, 2.0) - beta.get_value(args) / dx + gamma.get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = theta * dt * (beta.get_value(args) / dx - 2*alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (beta.get_value(args) / dx - 2*alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = theta * dt * (alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    // i = N

    double MatrixSystem::Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = -theta * dt * (alpha.get_value(args) / pow(dx, 2.0) + beta.get_value(args) / dx + gamma.get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::a_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (alpha.get_value(args) / pow(dx, 2.0) + beta.get_value(args) / dx + gamma.get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = -theta * dt * (beta.get_value(args) / dx + 2*alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = -(1 - theta) * dt * (beta.get_value(args) / dx + 2*alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = theta * dt * (alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        std::vector<double> args = {sigma, r};
        double res = (1 - theta) * dt * (alpha.get_value(args) / pow(dx, 2.0));
        return res;
    }
}
