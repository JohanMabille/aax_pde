#include "closed_form.hpp"
#include "payoff.h"
#include <vector>
#include <iostream>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"

int main(int argc, const char * argv[])
{
    // test of matrix
    std::vector<double> a = {1., 2., 3.};
    std::vector<double> b = {4., 5., 6.};

    std::vector<std::vector<double>> c = {a, b};

    Eigen::MatrixXd m(3, 2);

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            m(i, j) = c[j][i];
        }
    }

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            std::cout << m(i, j) << ",";
        }
        std::cout << std::endl;
    }
}
