#include "MatrixSystem.h"
#include "CoefEquation.h"
#include "BoundaryCondition.h"
#include <vector>
#include <cmath>
#include "eigen-3.4.0/Eigen/Dense"

namespace system_matrix
{
    //i != 0, i!=N

    double MatrixSystem::Omega(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (2*alpha.get_value({sigma, r}) / pow(dx, 2.0) - gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double MatrixSystem::a_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (-2*alpha.get_value({sigma, r}) / pow(dx, 2.0) + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double MatrixSystem::b_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) - beta.get_value({sigma, r}) / 2*dx);
        return res;
    }

    double MatrixSystem::c_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) + beta.get_value({sigma, r}) / 2*dx);
        return res;
    }

    double MatrixSystem::d_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) + beta.get_value({sigma, r}) / 2*dx);
        return res;
    }

    double MatrixSystem::e_i(double i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) - beta.get_value({sigma, r}) / 2*dx);
        return res;
    }

    //i = 0

    double MatrixSystem::Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) - beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double MatrixSystem::a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) - beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double MatrixSystem::b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (beta.get_value({sigma, r}) / dx - 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (beta.get_value({sigma, r}) / dx - 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    // i = N

    double MatrixSystem::Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) + beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double MatrixSystem::a_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) + beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double MatrixSystem::b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (beta.get_value({sigma, r}) / dx + 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -(1 - theta) * dt * (beta.get_value({sigma, r}) / dx + 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double MatrixSystem::e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    MatrixSystem::MatrixSystem(coef_eq::CoefEquation alpha,
                     coef_eq::CoefEquation beta,
                     coef_eq::CoefEquation gamma,
                     coef_eq::CoefEquation delta,
                     double theta,
                     double dt,
                     double dx,
                     double sigma,
                     double r,
                     boundary::BoundaryCondition boundary_small_spot,
                     boundary::BoundaryCondition boundary_big_spot,
                     double time,
                     std::vector<double> Xt1)
    {
        int N = Xt1.size();

        Eigen::MatrixXd A_prime(N,N);
        Eigen::MatrixXd A_second(N,N);
        Eigen::MatrixXd b(N,1);

        // ######### body #########
        // A'
        for (int i=1; i < N-1; ++i)
        {
            for (int j=0; j < i-1; ++j)
            {
                A_prime(i, j) = 0;
            }

            A_prime(i, i-1) = e_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma) / Omega(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_prime(i, i) = a_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma) / Omega(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_prime(i, i+1) = d_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma) / Omega(i, theta, dt, dx, sigma, r, alpha, beta, gamma);

            for (int j=i+2; j < N; ++j)
            {
                A_prime(i, j) = 0;
            }
        }

        // A''
        for (int i=1; i < N-1; ++i)
        {
            for (int j=0; j < i-1; ++j)
            {
                A_second(i, j) = 0;
            }

            A_second(i, i-1) = b_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma) / Omega(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_second(i, i) = 0;
            A_second(i, i+1) = c_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma) / Omega(i, theta, dt, dx, sigma, r, alpha, beta, gamma);

            for (int j=i+2; j < N; ++j)
            {
                A_second(i, j) = 0;
            }
        }

        // b'
        for (int i=0; i < N; ++i)
        {
            b(i,1) = delta.get_value({sigma, r});
        }

        // ######### boundaries #########
        // find a way to know if the condition is on f or the derivatives

        double small_spot = 0;
        double big_spot = 1000;
        std::vector<std::vector<double>> cond_small = boundary_small_spot.get_conditions(time, small_spot, N);
        std::vector<std::vector<double>> cond_big = boundary_big_spot.get_conditions(time, big_spot, N);
    }

    std::vector<double> MatrixSystem::get_result()
    {
        return m_result;
    }



}
