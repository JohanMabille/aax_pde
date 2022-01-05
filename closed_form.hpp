#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP

#include <vector>

namespace dauphine
{
    double vanilla_payoff(double fwd, double strike, bool is_call);
    double bs_time_value(double fwd, double strike, double volatility, double maturity);
    double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call);

    double call_delta(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double call_gamma(const double fwd, const double strike, const double r, const double volatility, const double maturity) ;
    double call_vega(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double call_theta(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double call_rho(const double fwd, const double strike, const double r, const double volatility, const double maturity);

    double put_delta(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double put_gamma(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double put_vega(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double put_theta(const double fwd, const double strike, const double r, const double volatility, const double maturity);
    double put_rho(const double fwd, const double strike, const double r, const double volatility, const double maturity);

    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity);
    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call);
}

#endif
