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


    MatrixSystem::MatrixSystem(coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma, coef_eq::CoefEquation delta,
                     double theta, double dt, double dx, double sigma, double r, double time,
                     boundary::BoundaryCondition boundary_small_spot, boundary::BoundaryCondition boundary_big_spot,
                     std::vector<double> Xt1)
    {
        int N = Xt1.size();

        Eigen::MatrixXd A_prime(N,N);
        Eigen::MatrixXd A_second(N,N);
        Eigen::MatrixXd Omega_matrix(N,N);
        Eigen::MatrixXd b(N,1);
        Eigen::VectorXd Xt1_matrix = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Xt1.data(), Xt1.size());

        // ######### boundaries #########

        double small_spot = 0;
        double big_spot = 1000;
        std::vector<std::vector<double>> cond_small = boundary_small_spot.get_conditions(time, small_spot, N, theta,  dt,  dx,  sigma,  r,
                                                        alpha, beta, gamma);
        std::vector<std::vector<double>> cond_big = boundary_big_spot.get_conditions(time, big_spot, N, theta,  dt,  dx,  sigma,  r,
                                                        alpha, beta, gamma);

        for (int j=0; j < N; ++j)
        {
            A_prime(0,j) = cond_small[0][j];
            A_prime(N,j) = cond_big[0][j];
            A_second(0,j) = cond_small[1][j];
            A_second(N,j) = cond_big[1][j];
        }

        Omega_matrix(0,0) = cond_small[2][0];
        Omega_matrix(N,N) = cond_big[2][0];

        // ######### body #########
        // Omega
        for (int i=1; i<N-1; ++i)
        {
            Omega_matrix(i,i) = Omega(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma);
        }

        // A'
        for (int i=1; i < N-1; ++i)
        {
            for (int j=0; j < i-1; ++j)
            {
                A_prime(i, j) = 0;
            }

            A_prime(i, i-1) = e_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_prime(i, i) = a_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_prime(i, i+1) = d_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma);

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

            A_second(i, i-1) = b_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_second(i, i) = 0;
            A_second(i, i+1) = c_i(static_cast<double>(i), theta, dt, dx, sigma, r, alpha, beta, gamma);

            for (int j=i+2; j < N; ++j)
            {
                A_second(i, j) = 0;
            }
        }

        // b'
        for (int i=0; i < N; ++i)
        {
            b(i,1) = delta.get_value({sigma, r}) * dt;
        }

        // ######### Computations to set up system #########
        m_A = Omega_matrix - A_second;
        m_b = A_prime * Xt1_matrix + b;
    }

    std::vector<double> MatrixSystem::solve()
    {
        Eigen::VectorXd res_vec = m_A.inverse() * m_b;

        std::vector<double> result;
        result.resize(res_vec.size());
        Eigen::VectorXd::Map(&result[0], res_vec.size()) = res_vec;

        return result;
    }


}
