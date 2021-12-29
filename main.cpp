#include "closed_form.hpp"
#include "payoff.h"
#include "Mesh.h"
#include <vector>
#include <iostream>
#include "eigen-3.4.0/Eigen/Dense"

void test_mesh()
{
    // params
    double S = 100.;
    double K = 100.;
    double sigma = 0.16;
    double theta = 0.5;
    int maturity = 30;
    int nb_steps_time = maturity;
    int nb_steps_space = 11;
    double r = 0.;
    payoff::Payoff *pf = new payoff::Call(K);
    boundary::BoundaryCondition *b_small = new boundary::ConditionSmall();
    boundary::BoundaryCondition *b_big = new boundary::ConditionBig();
    coef_eq::CoefEquation *alpha = new coef_eq::Alpha();
    coef_eq::CoefEquation *beta = new coef_eq::Beta();
    coef_eq::CoefEquation *gamma = new coef_eq::Gamma();
    coef_eq::CoefEquation *delta = new coef_eq::Delta();

    // pricing pde
    std::cout << "start" << std::endl;
    mesh::Mesh Mesh_call(S, sigma, maturity, nb_steps_space, nb_steps_time, theta, r, pf, b_small, b_big, *alpha, *beta, *gamma, *delta);
    Mesh_call.run();
    std::cout << Mesh_call.get_price() << std::endl;

    // pricing closed form


    delete pf;
    delete b_small;
    delete b_big;
}

int main(int argc, const char * argv[])
{
    test_mesh();
}
