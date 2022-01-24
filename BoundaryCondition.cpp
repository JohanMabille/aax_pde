#include "BoundaryCondition.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "payoff.h"

namespace boundary
{
    std::vector<std::vector<double>> BoundaryCondition::get_conditions(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                                        double spot_min, double spot_max, double S0, double maturity) const
    {
        std::vector<double> coef_Xt = get_coef_Xt(time, space, length, theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity);
        std::vector<double> coef_Xt1 = get_coef_Xt1(time, space, length, theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity);
        std::vector<double> coef_diag = get_coef_diag(time, space, length, theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity);
        std::vector<std::vector<double>> conditions = {coef_Xt, coef_Xt1, coef_diag};
        return conditions;
    }

    double ConditionBig::Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = -theta * dt * (alpha->get_value(args) / pow(dx, 2.0) + beta->get_value(args) / dx + gamma->get_value(args)) - 1;
        return res;
    }

    double ConditionBig::a_N(double theta, double
                              dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                               double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = (1 - theta) * dt * (alpha->get_value(args) / pow(dx, 2.0) + beta->get_value(args) / dx + gamma->get_value(args)) - 1;
        return res;
    }

    double ConditionBig::b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                              double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = -theta * dt * (beta->get_value(args) / dx + 2*alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionBig::c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                              double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = -(1 - theta) * dt * (beta->get_value(args) / dx + 2*alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionBig::d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                              double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = theta * dt * (alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionBig::e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                              double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = (1 - theta) * dt * (alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                    double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = -theta * dt * (alpha->get_value(args) / pow(dx, 2.0) - beta->get_value(args) / dx + gamma->get_value(args)) - 1;
        return res;
    }

    double ConditionSmall::a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = (1 - theta) * dt * (alpha->get_value(args) / pow(dx, 2.0) - beta->get_value(args) / dx + gamma->get_value(args)) - 1;
        return res;
    }

    double ConditionSmall::b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = theta * dt * (beta->get_value(args) / dx - 2*alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = (1 - theta) * dt * (beta->get_value(args) / dx - 2*alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = theta * dt * (alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                double spot_min, double spot_max, double S0, double maturity, int time, int space) const
    {
        std::vector<double> args = {sigma, r, dx, dt, maturity, S0, spot_min, spot_max, (double)time, (double)space};
        double res = (1 - theta) * dt * (alpha->get_value(args) / pow(dx, 2.0));
        return res;
    }

    std::vector<double> ConditionSmall::get_coef_Xt(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, double spot_min, double spot_max, double S0, double maturity) const
                                                                                                        // (zeros, 0, b0, d0)
    {
        std::vector<double> res(length, 0.0);
        res[length - 2] = ConditionSmall::b_0(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        res[length - 1] = ConditionSmall::d_0(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);

        return res;
    }


    std::vector<double> ConditionSmall::get_coef_Xt1(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, double spot_min, double spot_max, double S0, double maturity) const
                                                    // (zeros, a0, c0, e0)
    {
        std::vector<double> res(length, 0.0);
        res[length - 3] = ConditionSmall::a_0(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        res[length - 2] = ConditionSmall::c_0(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        res[length - 1] = ConditionSmall::e_0(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);

        return res;
    }


    std::vector<double> ConditionSmall::get_coef_diag(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, double spot_min, double spot_max, double S0, double maturity) const
                                                                                                            // (zeros,Omega0)
    {
        std::vector<double> res(length, 0.0);
        res[length - 1] = ConditionSmall::Omega_0(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        return res;
    }

    std::vector<double> ConditionBig::get_coef_Xt(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, double spot_min, double spot_max, double S0, double maturity) const
                                                    // (dN, bN, 0, zeros)
    {
        std::vector<double> res(length, 0.0);
        res[0] = ConditionBig::d_N(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        res[1] = ConditionBig::b_N(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        return res;
    }


    std::vector<double> ConditionBig::get_coef_Xt1(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, double spot_min, double spot_max, double S0, double maturity) const
                                                    // (eN, cN, aN, zeros)
    {
        std::vector<double> res(length, 0.0);
        res[0] = ConditionBig::e_N(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        res[1] = ConditionBig::c_N(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        res[2] = ConditionBig::a_N(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        return res;
    }


    std::vector<double> ConditionBig::get_coef_diag(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma, double spot_min, double spot_max, double S0, double maturity) const
                                                    // (OmegaN, zeros)
    {
        std::vector<double> res(length, 0.0);
        res[0] = ConditionBig::Omega_N(theta, dt, dx, sigma, r, alpha, beta, gamma, spot_min, spot_max, S0, maturity, time, space);
        return res;
    }
}
