#ifndef COEFEQUATION_H_INCLUDED
#define COEFEQUATION_H_INCLUDED
#include <vector>

namespace coef_eq
{
    // The idea of being able to overload coefficient is good, but
    // it coule be improved:
    // - all the coefficientd could be gathered in a single class
    // (PDEModel for instance)
    // - the api could be a virtual method per coefficient that fills
    // the coefficient for the whole grid:
    // void get_alpha(std::vector<std::vector<double>>& mesh, args...)
    // This way you avoid multiple call to virtual functions in the
    // matrix computation and improve performance: the access operator
    // of std::vector can be inlined, while a call to a virtual method
    // cannnot. This prevents a lot of optimizations a compiler can do.
    class CoefEquation
    {
    public:
        // Why not taking the argument by const ref to avoid a copy? (
        // and an underlying dynamic memory allocation / deletion)
        virtual double get_value(const std::vector<double>& args) const = 0;
        virtual ~CoefEquation() = default;
        // Almost good:
        CoefEquation(const CoefEquation&) = delete;
        CoefEquation& operator=(const CoefEquation&) = delete;
        CoefEquation(CoefEquation&&) = delete;
        CoefEquation& operator=(CoefEquation&&) = delete;

    protected:
        CoefEquation() = default;
    };

    class Alpha: public CoefEquation
    {
    public:
        Alpha() = default;
        double get_value(const std::vector<double>& args) const override; //Override: pour utiliser le même nom de méthode sur la fille
    };

    class Beta: public CoefEquation
    {
    public:
        Beta() = default;
        double get_value(const std::vector<double>& args) const override;
    };

    class Gamma: public CoefEquation
    {
    public:
        Gamma() = default;
        double get_value(const std::vector<double>& args) const override;
    };

    class Delta: public CoefEquation
    {
    public:
        Delta() = default;
        double get_value(const std::vector<double>& args) const override;
    };



}

#endif // COEFEQUATION_H_INCLUDED
