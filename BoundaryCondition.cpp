#include "BoundaryCondition.h"
#include <vector>
#include <iostream>
#include <cmath>

namespace boundary
{
    std::vector<std::vector<double>> BoundaryCondition::get_conditions(int time, double space, int length) const
    {
        std::vector<double> coef_Xt = BoundaryCondition::get_coef_Xt(time, space, length);
        std::vector<double> coef_Xt1 = BoundaryCondition::get_coef_Xt1(time, space, length);
        std::vector<double> coef_diag = BoundaryCondition::get_coef_diag(time, space, length);
        std::vector<std::vector<double>> conditions = {coef_Xt, coef_Xt1, coef_diag};
        return conditions;
    }

    std::vector<double> BoundaryCondition::get_coef_Xt(int time, double space, int length) const
    {
        //todo
        std::vector<double> condition;
        return condition;
    }

    std::vector<double> BoundaryCondition::get_coef_Xt1(int time, double space, int length) const
    {
        //todo
        std::vector<double> condition;
        return condition;
    }

    std::vector<double> BoundaryCondition::get_coef_diag(int time, double space, int length) const
    {
        //todo
        std::vector<double> condition;
        return condition;
    }

    //TODO: implementer un throw dans le cas où les variables ne sont pas déclarées

    double ConditionBig::Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) + beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double ConditionBig::a_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) + beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double ConditionBig::b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (beta.get_value({sigma, r}) / dx + 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionBig::c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -(1 - theta) * dt * (beta.get_value({sigma, r}) / dx + 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionBig::d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionBig::e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) - beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double ConditionSmall::a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0) - beta.get_value({sigma, r}) / dx + gamma.get_value({sigma, r})) - 1;
        return res;
    }

    double ConditionSmall::b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (beta.get_value({sigma, r}) / dx - 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (beta.get_value({sigma, r}) / dx - 2*alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }

    double ConditionSmall::e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha.get_value({sigma, r}) / pow(dx, 2.0));
        return res;
    }
}
