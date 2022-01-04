#ifndef COEFEQUATION_H_INCLUDED
#define COEFEQUATION_H_INCLUDED
#include <vector>

namespace coef_eq
{
    class CoefEquation
    {
    public:
        virtual double get_value(std::vector<double> args) const = 0;
    };

    class Alpha: public CoefEquation
    {
    public:
        double get_value(std::vector<double> args) const override; //Override: pour utiliser le même nom de méthode sur la fille
    };

    class Beta: public CoefEquation
    {
    public:
        double get_value(std::vector<double> args) const override;
    };

    class Gamma: public CoefEquation
    {
    public:
        double get_value(std::vector<double> args) const override;
    };

    class Delta: public CoefEquation
    {
    public:
        double get_value(std::vector<double> args) const override;
    };



}

#endif // COEFEQUATION_H_INCLUDED
