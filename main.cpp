#include "closed_form.hpp"
#include "payoff.h"
#include "Mesh.h"
#include <vector>
#include <iostream>
#include "payoff.h"
#include "eigen-3.4.0/Eigen/Dense"
#include <cmath>

double roundoff(double value)
{
  double pow_10 = pow(10.0, 4.0);
  return round(value * pow_10) / pow_10;
}

void print_grid(std::vector<std::vector<double>> grid)
{
    for (int i = 0; i < (int)grid.size(); i++)
    {
        payoff::print_vector(grid[i]);
    }
}

void test_mesh()
{
    // params
    double S = 100.;
    double K = 120.;
    double sigma = 0.16;
    double theta = 0.5;
<<<<<<< HEAD

    double maturity = 1;
    int nb_steps_time = maturity*365;
    int nb_steps_space = 100.;
    double r = 0.05;


=======
    double maturity = 0.25;
    int nb_steps_time = maturity * 1000;
    int nb_steps_space = 501.;
    double r = 0.03;
>>>>>>> 71455175705a76c8ebc7861a3fc49364041309b8
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

    std::vector<std::vector<double>> mesh_viz = Mesh_call.get_mesh();
    double price_pde = Mesh_call.get_price();
    double delta_pde = Mesh_call.get_delta();
    double gamma_pde = Mesh_call.get_gamma();
    double vega_pde = Mesh_call.get_vega();
    double theta_pde = Mesh_call.get_theta();

    // pricing closed form
    double price_cf = dauphine::bs_price(S*std::exp(r*maturity), K, sigma, maturity, true);
    double delta_cf = dauphine::call_delta(S,K,r,sigma,maturity);
    double gamma_cf = dauphine::call_gamma(S,K,r,sigma,maturity);
    double vega_cf = dauphine::call_vega(S,K,r,sigma,maturity);
    double theta_cf = dauphine::call_theta(S,K,r,sigma,maturity);

//    std::cout << "\n ------------------GRID------------------ \n" << std::endl;
//    print_grid(mesh_viz);
    std::cout << "\n          PDE     ||   Closed form \n" << std::endl;
    std::cout << "Price:     " << roundoff(price_pde) << "   ||" << "   " << roundoff(price_cf) << std::endl;
    std::cout << "Delta:     " << roundoff(delta_pde) << "  ||" << "   " << roundoff(delta_cf) << std::endl;
    std::cout << "Gamma:     " << roundoff(gamma_pde) << "  ||" << "   " << roundoff(gamma_cf) << std::endl;
    std::cout << "Vega:      " << roundoff(vega_pde) << "   ||" << "   " << roundoff(vega_cf) << std::endl;
    std::cout << "Theta:     " << roundoff(theta_pde) << " ||" << "   " << roundoff(theta_cf) << std::endl;

    delete pf;
    delete b_small;
    delete b_big;
}

int main(int argc, const char * argv[])
{
    test_mesh();
}
