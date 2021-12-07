#include "MatrixSystem.h"
#include "CoefEquation.h"
#include "BoundaryCondition.h"
#include <vector>
#include "eigen-3.4.0/Eigen/Dense"

namespace system_matrix
{
    MatrixSystem::MatrixSystem(coef_eq::CoefEquation alpha,
                     coef_eq::CoefEquation beta,
                     coef_eq::CoefEquation gamma,
                     coef_eq::CoefEquation delta,
                     double theta,
                     boundary::BoundaryCondition boundary_small_spot,
                     boundary::BoundaryCondition boundary_big_spot,
                     std::vector<double> Xt1)
    {
        //todo
    }

    std::vector<double> MatrixSystem::get_result()
    {
        return m_result;
    }
}
