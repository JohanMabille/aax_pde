#include "BoundaryCondition.h"
#include <vector>
#include <iostream>

namespace boundary
{
    std::vector<std::vector<double>> BoundaryCondition::get_conditions(int time, double space)
    {
        std::vector<double> coef_Xt = BoundaryCondition::get_coef_Xt(time, space);
        std::vector<double> coef_Xt1 = BoundaryCondition::get_coef_Xt1(time, space);
        std::vector<std::vector<double>> conditions = {coef_Xt, coef_Xt1};
        return conditions;
    }

    std::vector<double> BoundaryCondition::get_coef_Xt(int time, double space)
    {
        //todo
        std::vector<double> condition;
        return condition;
    }

    std::vector<double> BoundaryCondition::get_coef_Xt1(int time, double space)
    {
        //todo
        std::vector<double> condition;
        return condition;
    }
}
