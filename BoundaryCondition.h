#ifndef BOUNDARYCONDITION_H_INCLUDED
#define BOUNDARYCONDITION_H_INCLUDED
# include <vector>
#include <stdexcept>

namespace boundary
{

    class BoundaryCondition
    {
    public:
        BoundaryCondition() = default;
        ~BoundaryCondition() = default;
        std::vector<std::vector<double>> get_conditions(int time, double space, int length) const; //calcul des conditions aux bornes à chaque itération
        //Potential TODO: calculer les conditions aux bornes d'un coup

    protected:
        virtual std::vector<double> get_coef_Xt(int time, double space, int length) const;
        virtual std::vector<double> get_coef_Xt1(int time, double space, int length) const;
    };

    class Dirichlet_test: public BoundaryCondition
    {
    protected:
        virtual std::vector<double> get_coef_Xt(int time, double space, int length) const override;
        virtual std::vector<double> get_coef_Xt1(int time, double space, int length) const override;
    };
}


#endif // BOUNDARYCONDITION_H_INCLUDED
