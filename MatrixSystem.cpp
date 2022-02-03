#include "MatrixSystem.h"
#include "CoefEquation.h"
#include "BoundaryCondition.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "eigen-3.4.0/Eigen/Dense"

namespace system_matrix
{
    double MatrixSystem::Omega(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        // Gathering arguments of the same type but totally unrelated into a container to ivoke a function
        // is usually not recommended. First, it is hard to read and to maintain, second it is not optimal
        // since this results in a dynamic allocation / deletion. A structure for storing the arguments
        // would be more adapted
        std::vector<double> args = {sigma, r, dx, dt, m_maturity, m_S0, m_spot_min, m_spot_max, (double)m_loop, (double)i};
        double res = theta * dt * (2*alpha->get_value(args) / pow(dx, 2.0) - gamma->get_value(args)) - 1;
        return res;
    }


    double MatrixSystem::a_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        std::vector<double> args = {sigma, r, dx, dt, m_maturity, m_S0, m_spot_min, m_spot_max, (double)m_loop, (double)i};
        double res = (1 - theta) * dt * (-2*alpha->get_value(args) / pow(dx, 2.0) + gamma->get_value(args)) - 1;
        return res;
    }

    double MatrixSystem::b_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        std::vector<double> args = {sigma, r, dx, dt, m_maturity, m_S0, m_spot_min, m_spot_max, (double)m_loop, (double)i};
        double res = theta * dt * (alpha->get_value(args) / pow(dx, 2.0) - beta->get_value(args) / (2*dx));
        return res;
    }

    double MatrixSystem::c_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        std::vector<double> args = {sigma, r, dx, dt, m_maturity, m_S0, m_spot_min, m_spot_max, (double)m_loop, (double)i};
        double res = theta * dt * (alpha->get_value(args) / pow(dx, 2.0) + beta->get_value(args) / (2*dx));
        return res;
    }

    double MatrixSystem::d_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        std::vector<double> args = {sigma, r, dx, dt, m_maturity, m_S0, m_spot_min, m_spot_max, (double)m_loop, (double)i};
        double res = (1 - theta) * dt * (alpha->get_value(args) / pow(dx, 2.0) + beta->get_value(args) / (2*dx));
        return res;
    }

    double MatrixSystem::e_i(int i, double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        std::vector<double> args = {sigma, r, dx, dt, m_maturity, m_S0, m_spot_min, m_spot_max, (double)m_loop, (double)i};
        double res = (1 - theta) * dt * (alpha->get_value(args) / pow(dx, 2.0) - beta->get_value(args) / (2*dx));
        return res;
    }



    MatrixSystem::MatrixSystem(coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, coef_eq::CoefEquation* delta,
                     double theta, double dt, double dx, double sigma, double r,
                     boundary::BoundaryCondition *boundary_small_spot, boundary::BoundaryCondition *boundary_big_spot,
                     std::vector<double> Xt1, double spot_min, double spot_max, double S0, int loop, double maturity):
                         m_spot_min(spot_min), m_spot_max(spot_max), m_S0(S0), m_loop(loop), m_maturity(maturity)
    {
        int N = Xt1.size();
        Eigen::MatrixXd A_prime = Eigen::MatrixXd::Zero(N, N);
        Eigen::MatrixXd A_second = Eigen::MatrixXd::Zero(N, N);
        Eigen::MatrixXd Omega_matrix = Eigen::MatrixXd::Zero(N, N);
        Eigen::MatrixXd b(N, 1);
        Eigen::VectorXd Xt1_matrix = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Xt1.data(), Xt1.size());

        // ######### boundaries #########

        std::vector<std::vector<double>> cond_small = boundary_small_spot -> get_conditions(loop, N, N, theta,  dt,  dx,  sigma,  r,
                                                        alpha, beta, gamma, spot_min, spot_max, S0, maturity);
        std::vector<std::vector<double>> cond_big = boundary_big_spot -> get_conditions(loop, 0, N, theta,  dt,  dx,  sigma,  r,
                                                        alpha, beta, gamma, spot_min, spot_max, S0, maturity);

        for (int j=0; j < N; ++j)
        {
            A_second(0,j) = cond_big[0][j];
            A_second(N-1,j) = cond_small[0][j];
            A_prime(0,j) = cond_big[1][j];
            A_prime(N-1,j) = cond_small[1][j];

        }

        // ######### body #########
        // Omega

        for (int i=0; i<N; ++i)
        {
            Omega_matrix(i,i) = Omega(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
        }


        Omega_matrix(0,0) = cond_big[2][0];
        Omega_matrix(N-1,N-1) = cond_small[2][N-1];

        // A'
        for (int i=1; i < N-1; ++i)
        {
            A_prime(i, i-1) = e_i(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_prime(i, i) = a_i(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_prime(i, i+1) = d_i(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
        }

        // A''
        for (int i=1; i < N-1; ++i)
        {
            A_second(i, i-1) = b_i(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
            A_second(i, i+1) = c_i(i, theta, dt, dx, sigma, r, alpha, beta, gamma);
        }

        // b'

        for (int i=0; i < N; ++i)
        {
            std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)loop, (double)i};
            b(i,0) = delta->get_value(args) * dt;
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
