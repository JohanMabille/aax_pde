#include "Mesh.h"
#include <cmath>
#include <iostream>

namespace mesh
{
    Mesh::Mesh(double& S0, double& sigma, int& maturity, int& nb_steps_space, int& nb_steps_time, double& theta, double& r,
             payoff::Payoff*& pf, boundary::BoundaryCondition*& bound_small, boundary::BoundaryCondition*& bound_big,
             coef_eq::CoefEquation& alpha, coef_eq::CoefEquation& beta, coef_eq::CoefEquation& gamma, coef_eq::CoefEquation& delta)
             : m_S0(S0), m_sigma(sigma), m_maturity(maturity), m_nb_steps_space(nb_steps_space), m_nb_steps_time(nb_steps_time), m_theta(theta), m_r(r),
             m_pf(pf), m_bound_small(bound_small), m_bound_big(bound_big), m_alpha(alpha), m_beta(beta), m_gamma(gamma), m_delta(delta)
             {}

    std::vector<double> Mesh::initiate_spot_values(double S0, double sigma, int maturity, int nb_steps)
    {
        std::cout << "inside initiate spot values" << std::endl;


        double spot_max = log(S0) + 5 * sigma * sqrt(maturity);
        double spot_min = log(S0) - 5 * sigma * sqrt(maturity);
        m_dx = (spot_max - spot_min) / nb_steps;

        std::vector<double> res(nb_steps);
        for (int i=0; i<nb_steps; ++i)
        {
            res[i] = spot_min + i * m_dx;
        }

        return res;
    }

    void Mesh::run()
    {
        std::cout << "inside Mesh run" << std::endl;

        std::vector<double> spot_axis = initiate_spot_values(m_S0, m_sigma, m_maturity, m_nb_steps_space);
        std::cout << "initiate_spot_values done" << std::endl;

        std::vector<double> Xt1 = m_pf->compute_payoff(spot_axis);
        std::cout << "compute_payoff done" << std::endl;

        m_dt = m_maturity / m_nb_steps_time;

        std::cout << "start loop" << std::endl;
        for (int i=0; i<m_nb_steps_time; ++i)
        {
            std::cout << "     loop number " << i << std::endl;
            double time = 0.0; // value test TODO: harmonize with the potential value needed in BoundaryCondition
            system_matrix::MatrixSystem matrix_system(m_alpha, m_beta, m_gamma, m_delta, m_theta, m_dt, m_dx, m_sigma, m_r, time, m_bound_small, m_bound_big, Xt1);
            std::cout << "     matrix system initiated" << std::endl;
            Xt1 = matrix_system.solve();
            std::cout << "     matrix system solved" << std::endl;
            grid_res.push_back(Xt1);
            std::cout << "     result pushed" << std::endl;
        }

        std::cout << "end of Mesh run" << std::endl;

    }

    double Mesh::get_price()
    {
        std::cout << "inside get_price" << std::endl;
        if (m_nb_steps_space % 2 == 1)
        {
            return grid_res[m_nb_steps_time - 1][floor(m_nb_steps_space / 2)];
        }
        else
        {
            // approximation if parameters given by the user are not perfect
            return (grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2] + grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2 - 1]) / 2;
        }
    }

    double Mesh::get_delta()
    {
        if (m_nb_steps_space % 2 == 1)
        {
            return (grid_res[m_nb_steps_time - 1][floor(m_nb_steps_space / 2) + 1] - grid_res[m_nb_steps_time - 1][floor(m_nb_steps_space / 2) - 1]) / (2 * m_dx);
        }
        else
        {
            // approximation if parameters given by the user are not perfect
            return (grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2 ] - grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2 - 1]) / m_dx;
        }
    }

    double Mesh::get_gamma()
    {
        if (m_nb_steps_space % 2 == 1)
        {
            return (grid_res[m_nb_steps_time - 1][floor(m_nb_steps_space / 2) + 2] - grid_res[m_nb_steps_time - 1][floor(m_nb_steps_space / 2) - 2]) / pow(2. * m_dx, 2.);
        }
        else
        {
            // approximation if parameters given by the user are not perfect
            return (grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2 + 1] - grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2 - 2]) / pow(2. * m_dx, 2.);
        }
    }

    double Mesh::get_theta()
    {
        if (m_nb_steps_space % 2 == 1)
        {
            return (grid_res[m_nb_steps_time - 1][floor(m_nb_steps_space / 2)] - grid_res[m_nb_steps_time - 2][floor(m_nb_steps_space / 2)]) / m_dt;
        }
        else
        {
            // approximation if parameters given by the user are not perfect
            return ((grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2] - grid_res[m_nb_steps_time - 1][m_nb_steps_space / 2 - 1]) - (grid_res[m_nb_steps_time - 2][m_nb_steps_space / 2] + grid_res[m_nb_steps_time - 2][m_nb_steps_space / 2 - 1])) / m_dt;
        }
    }

    double Mesh::get_vega(double bumbed_sigma)
    {
        if (vega_computed)
        {
            // compute vega using already computed second pricing
            return 2.0;
        }
        else
        {
            // launch pricing with bumbed sigma and compute vega
            vega_computed = true;
            return 0;
        }
    }
}

