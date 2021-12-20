#include "BoundaryCondition.h"
#include <vector>
#include <iostream>

namespace boundary
{
    std::vector<std::vector<double>> BoundaryCondition::get_conditions(int time, double space, int length) const
    {
        std::vector<double> coef_Xt = BoundaryCondition::get_coef_Xt(time, space, length);
        std::vector<double> coef_Xt1 = BoundaryCondition::get_coef_Xt1(time, space, length);
        std::vector<std::vector<double>> conditions = {coef_Xt, coef_Xt1};
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

    //TODO: implementer un throw dans le cas où les variables ne sont pas déclarées

    std::vector<double> Dirichlet_test::get_coef_Xt(int time, double space, int length) const
    {
        // f / [coefs]
        std::vector<double> res(length + 1, 0.0);
        res[0] = space;
        return res;
    }

    std::vector<double> Dirichlet_test::get_coef_Xt1(int time, double space, int length) const
    {
        // f / [coefs]
        std::vector<double> res(length + 1, 0.0);
        res[0] = space;
        return res;
    }
}
