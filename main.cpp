#include "closed_form.hpp"
#include "payoff.h"
#include "Mesh.h"
#include <vector>
#include <iostream>
#include "payoff.h"
#include "eigen-3.4.0/Eigen/Dense"
#include <cmath>

void test_mesh()
{
    // params
    double S = 100.;
    double K = 80.;
    double sigma = 0.20;
    double theta = 0.5;
    double maturity = 1;
    int nb_steps_time = maturity*365;
    int nb_steps_space = 1001.;
    double r = 0.05;
    payoff::Payoff *pf = new payoff::Call(K);
    boundary::BoundaryCondition *b_small = new boundary::ConditionSmall();
    boundary::BoundaryCondition *b_big = new boundary::ConditionBig();
    coef_eq::CoefEquation *alpha = new coef_eq::Alpha();
    coef_eq::CoefEquation *beta = new coef_eq::Beta();
    coef_eq::CoefEquation *gamma = new coef_eq::Gamma();
    coef_eq::CoefEquation *delta = new coef_eq::Delta();

    // pricing pde
    mesh::Mesh Mesh_call(S, sigma, maturity, nb_steps_space, nb_steps_time, theta, r, pf, b_small, b_big, alpha, beta, gamma, delta);
    Mesh_call.run();
    double price_pde = Mesh_call.get_price();
    double delta_pde = Mesh_call.get_delta();
    double gamma_pde = Mesh_call.get_gamma();
    double vega_pde = Mesh_call.get_vega();
    double theta_pde = Mesh_call.get_theta();

    // pricing closed form
    double price_cf = dauphine::bs_price(S*std::exp(r*maturity), K, sigma, maturity, true);
    double delta_cf = 0.50; // to implement ?
    double gamma_cf = 0.1;
    double vega_cf = 3.2;
    double theta_cf = -3.1;

    std::cout << "          PDE     ||     Closed form" << std::endl;
    std::cout << "price     " << price_pde << "    ||" << "     " << price_cf << std::endl;
    std::cout << "delta     " << delta_pde << "    ||" << "     " << delta_cf << std::endl;
    std::cout << "gamma     " << gamma_pde << "    ||" << "     " << gamma_cf << std::endl;
    std::cout << "vega      " << vega_pde << "    ||" << "     " << vega_cf << std::endl;
    std::cout << "theta     " << theta_pde << "    ||" << "     " << theta_cf << std::endl;

    delete pf;
    delete b_small;
    delete b_big;
}

int main(int argc, const char * argv[])
{
    test_mesh();
}
