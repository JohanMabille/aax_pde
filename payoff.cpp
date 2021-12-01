#include "payoff.h"
#include <iostream>
#include <vector>
#include <algorithm>

namespace payoff
{
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
        double k = 100; //to be removed
        std::transform(S.begin(), S.end(), res.begin(), [k](double x) {return std::max(0., x-k);});
        return res;
    }



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
        double k = 100; //to be removed
        std::transform(S.begin(), S.end(), res.begin(), [k](double x) {return std::max(0., k-x);});
        return res;
    }
}
