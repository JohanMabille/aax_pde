#include "Mesh.h"
#include <cmath>
#include <iostream>
#include "payoff.h"

namespace mesh
{
    Mesh::Mesh(double& S0, double& sigma, double& maturity, int& nb_steps_space, int& nb_steps_time, double& theta, double& r,
             payoff::Payoff*& pf, boundary::BoundaryCondition*& bound_small, boundary::BoundaryCondition*& bound_big,
             coef_eq::CoefEquation*& alpha, coef_eq::CoefEquation*& beta, coef_eq::CoefEquation*& gamma, coef_eq::CoefEquation*& delta)
             : m_S0(S0), m_sigma(sigma), m_maturity(maturity), m_nb_steps_space(nb_steps_space), m_nb_steps_time(nb_steps_time), m_theta(theta), m_r(r),
             m_pf(pf), m_bound_small(bound_small), m_bound_big(bound_big), m_alpha(alpha), m_beta(beta), m_gamma(gamma), m_delta(delta)
             {}

    std::vector<double> Mesh::initiate_spot_values(double S0, double sigma, double maturity, int nb_steps)
    {
        double spot_max = log(S0) + 5 * sigma * sqrt(maturity);
        double spot_min = log(S0) - 5 * sigma * sqrt(maturity);

        m_dx = (spot_max - spot_min) / nb_steps;

        std::vector<double> res(nb_steps);
        for (int i=0; i<nb_steps; ++i)
        {
            res[i] = spot_max - i * m_dx;
        }
        return res;
    }

    void Mesh::run(bool bumped)
    {
        if (bumped)
        {
            grid_res = main_run();
        }
        else
        {
            grid_res_bumped_sigma = main_run();
        }
    }

    std::vector<std::vector<double>> Mesh::main_run()
    {
        std::vector<double> log_spot_axis = initiate_spot_values(m_S0, m_sigma, m_maturity, m_nb_steps_space);
        std::vector<double> spot_axis;
        spot_axis.reserve(log_spot_axis.size());
        for(size_t i = 0; i < log_spot_axis.size(); ++i)
        {
            spot_axis.push_back(std::exp(log_spot_axis[i]));
        }

        std::vector<double> Xt1 = m_pf->compute_payoff(spot_axis);

        m_dt = m_maturity / m_nb_steps_time;

        // grid_res.push_back(Xt1); //à rajouter comme la première colonne doit etre le payoff? inutilefffffffffffff

        for (int i=0; i<m_nb_steps_time; ++i)
        {
            std::cout << "Loop number: " << i << std::endl;
            double time = 0.0; // value test TODO: harmonize with the potential value needed in BoundaryCondition
            system_matrix::MatrixSystem matrix_system(m_alpha, m_beta, m_gamma, m_delta, m_theta, m_dt, m_dx, m_sigma, m_r, time, m_bound_small,
                                                      m_bound_big, Xt1, spot_axis[0], spot_axis[spot_axis.size() - 1]);

            Xt1 = matrix_system.solve();
            grid_res.push_back(Xt1);

        }

    }

    double Mesh::get_price(bool bumped)
    {
        if (bumped){ return get_price(grid_res_bumped_sigma);}
        else {return get_price(grid_res);}
    }

    double Mesh::get_price(std::vector<std::vector<double>> grid)
    {
        if (m_nb_steps_space % 2 == 1)
            {
                return grid[m_nb_steps_time - 1][floor(m_nb_steps_space / 2)];
            }
            else
            {
                // interpolation if parameters given by the user are not perfect
                double denom = std::exp(m_dx) - std::exp(-m_dx);
                double w1 = (std::exp(m_dx) - 1) / denom;
                double w2 = (1 - std::exp(-m_dx)) / denom;
                return grid[m_nb_steps_time - 1][m_nb_steps_space / 2 - 1] * w1 + grid[m_nb_steps_time - 1][m_nb_steps_space / 2] * w2;
            }
    }

    double Mesh::get_delta()
    {
        double denom = m_S0 * (std::exp(m_dx) - std::exp(-m_dx));
        double num;
        if (m_nb_steps_space % 2 == 1)
        {
            int index = floor(m_nb_steps_space / 2);
            num = grid_res[m_nb_steps_time - 1][index + 1] - grid_res[m_nb_steps_time - 1][index - 1];
        }
        else
        {
            int index = m_nb_steps_space / 2;
            num = grid_res[m_nb_steps_time - 1][index] - grid_res[m_nb_steps_time - 1][index - 1];
        }
        return num / denom;
    }

    double Mesh::get_gamma()
    {
        double h1 = m_S0 * std::exp(-m_dx);
        double h2 = m_S0 * std::exp(m_dx);
        double denom = h1 * h2 * (h2 - h1);
        double num;
        if (m_nb_steps_space % 2 == 1)
        {
            int index = floor(m_nb_steps_space / 2);
            num = h1 * grid_res[m_nb_steps_time - 1][index + 1] - h1 * h2 * grid_res[m_nb_steps_time - 1][index] + h2 * grid_res[m_nb_steps_time - 1][index - 1];
        }
        else
        {
            int index = m_nb_steps_space / 2;
            num = h1 * grid_res[m_nb_steps_time - 1][index] - h1 * h2 * get_price() + h2 * grid_res[m_nb_steps_time - 1][index - 1];
        }
        return num / denom;
    }

    double Mesh::get_theta()
    {
        double scale_to_daily = m_nb_steps_time / 365;
        if (m_nb_steps_space % 2 == 1)
        {
            return (grid_res[m_nb_steps_time - 2][floor(m_nb_steps_space / 2)] - get_price()) / m_dt * scale_to_daily;
        }
        else
        {
            // approximation if parameters given by the user are not perfect
            int index = m_nb_steps_space / 2;
            double denom = std::exp(m_dx) - std::exp(-m_dx);
            double w1 = (std::exp(m_dx) - 1) / denom;
            double w2 = (1 - std::exp(-m_dx)) / denom;
            double price_t1 = grid_res[m_nb_steps_time - 2][index - 1] * w1 + grid_res[m_nb_steps_time - 2][index] * w2;
            return (price_t1 - get_price()) / m_dt * scale_to_daily;
        }
    }

    double Mesh::get_vega()
    {
        if (!vega_computed)
        {
            // relaunch pricing with bumped vol
            run(true);
            vega_computed = true;
        }
        return get_price(true) - get_price() / 0.01; // 1% vol move

    }
}

