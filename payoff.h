#ifndef PAYOFF_H_INCLUDED
#define PAYOFF_H_INCLUDED
#include <vector>


namespace payoff
{

    class Payoff
    {
    public :
        virtual ~Payoff() = default;
        Payoff& operator=(const Payoff&) = delete;
        Payoff(const Payoff&&) = delete;
        Payoff& operator=(Payoff&) = delete;
        virtual std::vector<double> compute_payoff (std::vector<double>& spot) const = 0;


    protected:
        Payoff() = default;

    };

    class Call: public Payoff
    {
    public:
        Call(double K);
        ~Call() = default;
        std::vector<double> compute_payoff (std::vector<double>& spot) const override;

    private :
        double m_strike;
    };


    class Put: public Payoff
    {
    public:
        Put(double K);
        ~Put() = default;
        std::vector<double> compute_payoff (std::vector<double>& spot) const override;
    private :
        double m_strike;
    };

    void test_payoff();
    void print_vector(const std::vector<double>& v);


}
#endif // PAYOFF_H_INCLUDED
