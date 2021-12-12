#ifndef COEFEQUATION_H_INCLUDED
#define COEFEQUATION_H_INCLUDED

namespace coef_eq
{
    class CoefEquation
    {
    public:
        double alpha(double r, double sigma) const;
        double beta(double r, double sigma) const;
        double gamma(double r, double sigma) const;
        double delta(double r, double sigma) const;

        double Omega(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double a_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double b_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double c_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double d_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double e_i(double i, double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;

        double Omega_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double a_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double b_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double c_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double d_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double e_0(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;

        double Omega_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double a_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double b_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double c_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double d_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;
        double e_N(double theta, double dt, double dx, coef_eq::CoefEquation alpha, coef_eq::CoefEquation beta, coef_eq::CoefEquation gamma) const;




    };




}

#endif // COEFEQUATION_H_INCLUDED
