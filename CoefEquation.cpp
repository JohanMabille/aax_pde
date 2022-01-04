#include "CoefEquation.h"
#include <iostream>
#include <vector>
#include <cmath>

namespace coef_eq
{

    double Alpha::get_value(std::vector<double> args) const
    {
        double sigma = args[0];
        return -0.5*pow(sigma,2.0);
    }

    double Beta::get_value(std::vector<double> args) const
    {
        double sigma = args[0];
        double r = args[1];
        return 0.5*pow(sigma,2.0) + r;
    }

    double Gamma::get_value(std::vector<double> args) const
    {
        double r = args[1];
        return -r;
    }

    double Delta::get_value(std::vector<double> args) const
    {
        return 0.0;
    }

}

