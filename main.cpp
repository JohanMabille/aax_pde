#include "closed_form.hpp"
#include <iostream>
#include "payoff.h"
#include <iterator>


void print_vector(const std::vector<double>& v)
{
    std::ostream_iterator<double> out_it(std::cout, ", ");
    std::copy(v.begin(), v.end(), out_it);
    std::cout << std::endl;
}


int main(int argc, const char * argv[])
{
    std::vector<double> S = {98., 100., 102.};

    payoff::Call test_call(100);
    print_vector(test_call.compute_payoff(S));

    payoff::Put test_put(100);
    print_vector(test_put.compute_payoff(S));

}
