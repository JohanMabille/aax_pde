#ifndef BOUNDARYCONDITION_H_INCLUDED
#define BOUNDARYCONDITION_H_INCLUDED
# include <vector>
#include <stdexcept>
#include "CoefEquation.h"

namespace boundary
{

    class BoundaryCondition
    {
    public:
        BoundaryCondition() = default;
//        virtual ~BoundaryCondition() = 0;
        ~BoundaryCondition();
        std::vector<std::vector<double>> get_conditions(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const; //calcul des conditions aux bornes à chaque itération
        //Potential TODO: calculer les conditions aux bornes d'un coup

    protected:
        virtual std::vector<double> get_coef_Xt(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const =0;

        virtual std::vector<double> get_coef_Xt1(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const = 0;

        virtual std::vector<double> get_coef_diag(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const =0;
    };


    class ConditionSmall: public BoundaryCondition
    {

    public:
        ~ConditionSmall();

    protected:
        virtual std::vector<double> get_coef_Xt(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const override;

        virtual std::vector<double> get_coef_Xt1(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const override;

        virtual std::vector<double> get_coef_diag(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const override;

    private:
        double Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
    };


    class ConditionBig: public BoundaryCondition
    {

    public:
        ~ConditionBig();

    protected:
        virtual std::vector<double> get_coef_Xt(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const override;

        virtual std::vector<double> get_coef_Xt1(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const override;

        virtual std::vector<double> get_coef_diag(int time, double space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const override;
    private:
        double Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double a_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;
        double e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma) const;

    };
}


#endif // BOUNDARYCONDITION_H_INCLUDED
