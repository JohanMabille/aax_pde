#include "BoundaryCondition.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "payoff.h"

namespace boundary
{
    std::vector<std::vector<double>> BoundaryCondition::get_conditions(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        std::vector<double> coef_Xt = get_coef_Xt(time, space, length, theta, dt, dx, sigma, r, alpha, beta, gamma);
        std::vector<double> coef_Xt1 = get_coef_Xt1(time, space, length, theta, dt, dx, sigma, r, alpha, beta, gamma);
        std::vector<double> coef_diag = get_coef_diag(time, space, length, theta, dt, dx, sigma, r, alpha, beta, gamma);
        std::vector<std::vector<double>> conditions = {coef_Xt, coef_Xt1, coef_diag};
        return conditions;
    }

    BoundaryCondition::~BoundaryCondition()
    {
        std::cout << "BC destructor" << std::endl;
    }

    double ConditionBig::Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = -theta * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0) + beta->get_value({sigma, r}) / dx + gamma->get_value({sigma, r})) - 1;
        // std::cout <<"Omega_N: "<<res<<std::endl;
        return res;
    }

    double ConditionBig::a_N(double theta, double
                              dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = (1 - theta) * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0) + beta->get_value({sigma, r}) / dx + gamma->get_value({sigma, r})) - 1;
//        std::cout <<"aN: "<<res<<std::endl;
        return res;
    }

    double ConditionBig::b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = -theta * dt * (beta->get_value({sigma, r}) / dx + 2*alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"bN: "<<res<<std::endl;
        return res;
    }

    double ConditionBig::c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = -(1 - theta) * dt * (beta->get_value({sigma, r}) / dx + 2*alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"cN: "<<res<<std::endl;
        return res;
    }

    double ConditionBig::d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = theta * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"dN: "<<res<<std::endl;
        return res;
    }

    double ConditionBig::e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = (1 - theta) * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"eN: "<<res<<std::endl;
        return res;
    }

    double ConditionSmall::Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = -theta * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0) - beta->get_value({sigma, r}) / dx + gamma->get_value({sigma, r})) - 1;
//        std::cout <<"Omega_0: "<<res<<std::endl;
        return res;
    }

    double ConditionSmall::a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = (1 - theta) * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0) - beta->get_value({sigma, r}) / dx + gamma->get_value({sigma, r})) - 1;
//        std::cout <<"a_0: "<<res<<std::endl;
        return res;
    }

    double ConditionSmall::b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = theta * dt * (beta->get_value({sigma, r}) / dx - 2*alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"b_0: "<<res<<std::endl;
        return res;
    }

    double ConditionSmall::c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = (1 - theta) * dt * (beta->get_value({sigma, r}) / dx - 2*alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"c_0: "<<res<<std::endl;
        return res;
    }

    double ConditionSmall::d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = theta * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"d_0: "<<res<<std::endl;
        return res;
    }

    double ConditionSmall::e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
    {
        double res = (1 - theta) * dt * (alpha->get_value({sigma, r}) / pow(dx, 2.0));
//        std::cout <<"e_0: "<<res<<std::endl;
        return res;
    }



    ConditionSmall::~ConditionSmall()
    {
        std::cout << "BC Small destructor" << std::endl;
    }


    std::vector<double> ConditionSmall::get_coef_Xt(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
                                                                                                        // (zeros, 0, b0, d0)
    {
        std::vector<double> res(length, 0.0);
        res[length - 2] = ConditionSmall::b_0(theta, dt, dx, sigma, r, alpha, beta, gamma);
        res[length - 1] = ConditionSmall::d_0(theta, dt, dx, sigma, r, alpha, beta, gamma);

        return res;
    }


    std::vector<double> ConditionSmall::get_coef_Xt1(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const

                                                    // (zeros, a0, c0, e0)
    {
        std::vector<double> res(length, 0.0);
        res[length - 3] = ConditionSmall::a_0(theta, dt, dx, sigma, r, alpha, beta, gamma);
        res[length - 2] = ConditionSmall::c_0(theta, dt, dx, sigma, r, alpha, beta, gamma);
        res[length - 1] = ConditionSmall::e_0(theta, dt, dx, sigma, r, alpha, beta, gamma);

        return res;
    }

    // checker qu'on a bien mis les bons coef pour Xt et Xt1

    std::vector<double> ConditionSmall::get_coef_diag(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
                                                                                                            // (zeros,Omega0)
    {
        std::vector<double> res(length, 0.0);
        res[length - 1] = ConditionSmall::Omega_0(theta, dt, dx, sigma, r, alpha, beta, gamma);
        return res;
    }



    ConditionBig::~ConditionBig()
    {
        std::cout << "BC Big destructor" << std::endl;
    }

    std::vector<double> ConditionBig::get_coef_Xt(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
                                                    // (dN, bN, 0, zeros)
    {
        std::vector<double> res(length, 0.0);
        res[0] = ConditionBig::d_N(theta, dt, dx, sigma, r, alpha, beta, gamma);
        res[1] = ConditionBig::b_N(theta, dt, dx, sigma, r, alpha, beta, gamma);
        return res;
    }


    std::vector<double> ConditionBig::get_coef_Xt1(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
                                                    // (eN, cN, aN, zeros)
    {
        std::vector<double> res(length, 0.0);
        res[0] = ConditionBig::e_N(theta, dt, dx, sigma, r, alpha, beta, gamma);
        res[1] = ConditionBig::c_N(theta, dt, dx, sigma, r, alpha, beta, gamma);
        res[2] = ConditionBig::a_N(theta, dt, dx, sigma, r, alpha, beta, gamma);
        return res;
    }


    std::vector<double> ConditionBig::get_coef_diag(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const
                                                    // (OmegaN, zeros)
    {
        std::vector<double> res(length, 0.0);
        res[0] = ConditionBig::Omega_N(theta, dt, dx, sigma, r, alpha, beta, gamma);
        return res;
    }
}
