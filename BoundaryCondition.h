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
        virtual ~BoundaryCondition() = default;
        //Boundary conditions computations at each iteration
        std::vector<std::vector<double>> get_conditions(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const;

        BoundaryCondition& operator=(const BoundaryCondition&) = delete;
        BoundaryCondition(const BoundaryCondition&&) = delete;
        BoundaryCondition& operator=(BoundaryCondition&) = delete;

    protected:
        BoundaryCondition() = default;
        virtual std::vector<double> get_coef_Xt(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                double spot_min, double spot_max, double S0, double maturity) const =0;

        virtual std::vector<double> get_coef_Xt1(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                double spot_min, double spot_max, double S0, double maturity) const = 0;

        virtual std::vector<double> get_coef_diag(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                    coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                    double spot_min, double spot_max, double S0, double maturity) const =0;
    };


    class ConditionSmall: public BoundaryCondition
    {

    public:
        ConditionSmall() = default;

    protected:
        virtual std::vector<double> get_coef_Xt(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const override;

        virtual std::vector<double> get_coef_Xt1(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const override;

        virtual std::vector<double> get_coef_diag(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const override;

    private:
        double Omega_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                        double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double a_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                        double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double b_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                        double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double c_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                        double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double d_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                        double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double e_0(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                        double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
    };


    class ConditionBig: public BoundaryCondition
    {

    public:
        ConditionBig() = default;

    protected:
        virtual std::vector<double> get_coef_Xt(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const override;

        virtual std::vector<double> get_coef_Xt1(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const override;

        virtual std::vector<double> get_coef_diag(int time, int space, int length, double theta, double dt, double dx, double sigma, double r,
                                                        coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                                        double spot_min, double spot_max, double S0, double maturity) const override;
    private:
        double Omega_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double a_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double b_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double c_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double d_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const;
        double e_N(double theta, double dt, double dx, double sigma, double r, coef_eq::CoefEquation* alpha, coef_eq::CoefEquation* beta, coef_eq::CoefEquation* gamma,
                                 double spot_min, double spot_max, double S0, double maturity, int time, int space) const;

    };
}


#endif // BOUNDARYCONDITION_H_INCLUDED
