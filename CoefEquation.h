#ifndef COEFEQUATION_H_INCLUDED
#define COEFEQUATION_H_INCLUDED
#include <vector>

namespace coef_eq
{
    class CoefEquation
    {
    public:
        virtual double get_value(std::vector<double> args) const = 0;
        virtual ~CoefEquation() = default;
        CoefEquation& operator=(const CoefEquation&) = delete;
        CoefEquation(const CoefEquation&&) = delete;
        CoefEquation& operator=(CoefEquation&) = delete;

    protected:
        CoefEquation() = default;
    };

    class Alpha: public CoefEquation
    {
    public:
        Alpha() = default;
        double get_value(std::vector<double> args) const override; //Override: pour utiliser le même nom de méthode sur la fille
    };

    class Beta: public CoefEquation
    {
    public:
        Beta() = default;
        double get_value(std::vector<double> args) const override;
    };

    class Gamma: public CoefEquation
    {
    public:
        Gamma() = default;
        double get_value(std::vector<double> args) const override;
    };

    class Delta: public CoefEquation
    {
    public:
        Delta() = default;
        double get_value(std::vector<double> args) const override;
    };



}

#endif // COEFEQUATION_H_INCLUDED
