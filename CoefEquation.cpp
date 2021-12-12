#include "CoefEquation.h"
#include <iostream>

//i != 0, i!=N

namespace coef_eq
{


    double CoefEquation::Omega(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (2*alpha / (dx)^2 - gamma) - 1
        return res
    }

    double CoefEquation::a_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (-2*alpha / (dx)^2 + gamma) - 1
        return res
    }

    double CoefEquation::b_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha / (dx)^2 - beta / 2*dx)
        return res
    }

    double CoefEquation::c_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha / (dx)^2 + beta / 2*dx)
        return res
    }

    double CoefEquation::d_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha / (dx)^2 + beta / 2*dx)
        return res
    }

    double CoefEquation::e_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha / (dx)^2 - beta / 2*dx)
        return res
    }

    //i = 0

    double CoefEquation::Omega_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (alpha / (dx)^2 - beta / dx + gamma) - 1
        return res
    }

    double CoefEquation::a_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha / (dx)^2 - beta / dx - gamma) - 1
        return res
    }

    double CoefEquation::b_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (beta / dx - 2*alpha / (dx)^2)
        return res
    }

    double CoefEquation::c_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (beta / dx - 2*alpha / (dx)^2)
        return res
    }

    double CoefEquation::d_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha / (dx)^2)
        return res
    }

    double CoefEquation::e_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha / (dx)^2)
        return res
    }

    // i = N

    double CoefEquation::Omega_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (alpha / (dx)^2 + beta / dx + gamma) - 1
        return res
    }

    double CoefEquation::a_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha / (dx)^2 + beta / dx + gamma) - 1
        return res
    }

    double CoefEquation::b_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -theta * dt * (beta / dx + 2*alpha / (dx)^2)
        return res
    }

    double CoefEquation::c_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = -(1 - theta) * dt * (beta / dx + 2*alpha / (dx)^2)
        return res
    }

    double CoefEquation::d_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = theta * dt * (alpha / (dx)^2)
        return res
    }

    double CoefEquation::e_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const
    {
        double res = (1 - theta) * dt * (alpha / (dx)^2)
        return res
    }

}

