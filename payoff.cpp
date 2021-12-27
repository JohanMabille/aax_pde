#include "payoff.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>


namespace payoff
{

    // Call

    Call::Call(double K):m_strike(K)
    {
        std::cout << "Call constructor" << std::endl;
    }

    Call::~Call()
    {
        std::cout << "Call destructor" << std::endl;
    }


    std::vector<double> Call::compute_payoff(std::vector<double>& S) const
    {
        std::vector<double> res(S.size());
        std::transform(S.begin(), S.end(), res.begin(), [&](double x) {return std::max(0., x - m_strike);});
        return res;
    }

    // Put

    Put::Put(double K):m_strike(K)
    {
        std::cout << "Put constructor" << std::endl;
    }

    Put::~Put()
    {
        std::cout << "Put destructor" << std::endl;
    }


    std::vector<double> Put::compute_payoff(std::vector<double>& S) const
    {
        std::vector<double> res(S.size());
        std::transform(S.begin(), S.end(), res.begin(), [&](double x) {return std::max(0., m_strike - x);});
        return res;
    }


    // test

    void print_vector(const std::vector<double>& v)
    {
        std::ostream_iterator<double> out_it(std::cout, ", ");
        std::copy(v.begin(), v.end(), out_it);
        std::cout << std::endl;
    }

    void test_payoff()
    {
        std::vector<double> S = { 98., 100., 102. };

        payoff::Call test_call(100);
        print_vector(test_call.compute_payoff(S));

        payoff::Put test_put(100);
        print_vector(test_put.compute_payoff(S));
    }
}
