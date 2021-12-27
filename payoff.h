#ifndef PAYOFF_H_INCLUDED
#define PAYOFF_H_INCLUDED
#include <vector>


namespace payoff
{

    class Payoff
    {
    public :
        Payoff() = default;
        virtual std::vector<double> compute_payoff (std::vector<double>& spot) const = 0;
        ~Payoff() = default;
    };

    class Call: public Payoff
    {
    public:
        Call(double K);
        ~Call();
        std::vector<double> compute_payoff (std::vector<double>& spot) const override;
    private :
        double m_strike;
    };


    class Put: public Payoff
    {
    public:
        Put(double K);
        ~Put();
        std::vector<double> compute_payoff (std::vector<double>& spot) const override;
    private :
        double m_strike;
    };

    void test_payoff();

}
#endif // PAYOFF_H_INCLUDED
