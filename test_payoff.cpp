/*
#include <iostream>
#include "payoff.h"
#include <iterator>



void print_vector(const std::vector<double>& v)
{
    std::ostream_iterator<double> out_it(std::cout, ", ");
    std::copy(v.begin(), v.end(), out_it);
    std::cout << std::endl;
}

int test_payoff()
{
    std::vector<double> S = {98., 100., 102.};

    payoff::Call test_call(100);
    return test_call.compute_payoff(S);


}

*/
