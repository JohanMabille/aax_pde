#ifndef MATRIXSYSTEM_H_INCLUDED
#define MATRIXSYSTEM_H_INCLUDED

#include <Eigen/Dense>
#include <vector>
#include "BoundaryCondition.h"
#include "CoefEquation.h"

namespace system
{
    class MatrixSystem
    {
    public:
        MatrixSystem(coef_eq::CoefEquation alpha,
                     coef_eq::CoefEquation beta,
                     coef_eq::CoefEquation gamma,
                     coef_eq::CoefEquation delta,
                     boundary::BoundaryCondition boundary_small_spot,
                     boundary::BoundaryCondition boundary_big_spot,
                     std::vector<double> Xt1);
        ~MatrixSystem() = default;
        void solve();
        std::vector<double> get_result();

    private:
        eigen::matrixXd m_A_prime;
        eigen::matrixXd m_A_second;
        eigen::matrixXd m_b;
        std::vector<double> result; // use of vector in order to avoid having to import eigen in client
    };
}


#endif // MATRIXSYSTEM_H_INCLUDED
