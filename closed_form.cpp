#include "closed_form.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <iostream>



namespace dauphine
{

    double norm_cdf(double x)
    {
        return 0.5 * std::erfc(-x / std::sqrt(2));
    }

    double norm_pdf(const double x) {
      return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
    }



    double vanilla_payoff(double fwd, double strike, bool is_call)
    {
        return std::max(is_call ? fwd - strike : strike - fwd, 0.);
    }

    double bs_time_value(double fwd, double strike, double volatility, double maturity)
    {
        if(strike == 0.)
        {
            return 0.;
        }
        else
        {
            double stddev = volatility * std::sqrt(maturity);
            if(stddev == 0.)
            {
                return 0.;
            }
            double tmp = std::log(fwd / strike) / stddev;
            double d1 = tmp + 0.5 * stddev;
            double d2 = tmp - 0.5 * stddev;
            double res;
            if(fwd > strike)
            {
                res = strike * norm_cdf(-d2) - fwd * norm_cdf(-d1);
            }
            else
            {
                res = fwd * norm_cdf(d1) - strike * norm_cdf(d2);
            }
            if(res <= std::numeric_limits<double>::min())
            {
                res = 0.;
            }
            return res;
        }
    }

    double bs_price(double fwd, double strike, double volatility, double maturity, bool is_call)
    {
        return vanilla_payoff(fwd, strike, is_call) + bs_time_value(fwd, strike, volatility, maturity);
    }

    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call)
    {
        std::vector<double> res(fwd.size());
        for(std::size_t i = 0; i < fwd.size(); ++i)
        {
            res[i] = vanilla_payoff(fwd[i], strike, is_call);
        }
        return res;
    }

    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity)
    {
        std::vector<double> res(fwd.size());
        for(std::size_t i = 0; i < fwd.size(); ++i)
        {
            res[i] = bs_time_value(fwd[i], strike, volatility, maturity);
        }
        return res;
    }

    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, bool is_call)
    {
        std::vector<double> res(fwd.size());
        for(std::size_t i = 0; i < fwd.size(); ++i)
        {
            res[i] = bs_price(fwd[i], strike, volatility, maturity, is_call);
        }
        return res;
    }




    // Calculation of D1 and D2
    double d_j(const int j, const double S, const double K, const double r, const double v, const double T) {
      return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
    }


    // GREEKS CALL

    double call_delta(const double S, const double K, const double r, const double v, const double T) {
      return norm_cdf(d_j(1, S, K, r, v, T));
    }

    double call_gamma(const double S, const double K, const double r, const double v, const double T) {
      return norm_pdf(d_j(1, S, K, r, v, T))/(S*v*sqrt(T));
    }

    double call_vega(const double S, const double K, const double r, const double v, const double T) {
      return S*norm_pdf(d_j(1, S, K, r, v, T))*sqrt(T);
    }

    double call_theta(const double S, const double K, const double r, const double v, const double T) {
      return (-(S*norm_pdf(d_j(1, S, K, r, v, T))*v)/(2*sqrt(T))
        - r*K*exp(-r*T)*norm_cdf(d_j(2, S, K, r, v, T)))/365;
    }

    double call_rho(const double S, const double K, const double r, const double v, const double T) {
      return K*T*exp(-r*T)*norm_cdf(d_j(2, S, K, r, v, T));
    }


    // GREEKS PUT


    double put_delta(const double S, const double K, const double r, const double v, const double T) {
      return norm_cdf(d_j(1, S, K, r, v, T)) - 1;
    }

    double put_gamma(const double S, const double K, const double r, const double v, const double T) {
      return call_gamma(S, K, r, v, T); // Identical to call by put-call parity
    }

    double put_vega(const double S, const double K, const double r, const double v, const double T) {
      return call_vega(S, K, r, v, T); // Identical to call by put-call parity
    }

    double put_theta(const double S, const double K, const double r, const double v, const double T) {
      return (-(S*norm_pdf(d_j(1, S, K, r, v, T))*v)/(2*sqrt(T))
        + r*K*exp(-r*T)*norm_cdf(-d_j(2, S, K, r, v, T)))/365;
    }

    double put_rho(const double S, const double K, const double r, const double v, const double T) {
      return -T*K*exp(-r*T)*norm_cdf(-d_j(2, S, K, r, v, T));
    }

}

