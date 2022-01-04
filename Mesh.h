#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include <vector>
#include "BoundaryCondition.h"
#include "MatrixSystem.h"
#include "payoff.h"
#include "CoefEquation.h"


namespace mesh
{
    class Mesh
    {

    public:
        Mesh(double& S0, double& sigma, double& maturity, int& nb_steps_space, int& nb_steps_time, double& theta, double& r,
             payoff::Payoff*& pf, boundary::BoundaryCondition*& bound_small, boundary::BoundaryCondition*& bound_big,
             coef_eq::CoefEquation*& alpha, coef_eq::CoefEquation*& beta, coef_eq::CoefEquation*& gamma, coef_eq::CoefEquation*& delta);
        ~Mesh () = default;
        void run(bool bumped = false);
        double get_price(bool bumped = false);
        double get_delta();
        double get_gamma();
        double get_theta();
        double get_vega();

    private:
        double m_S0;
        double m_sigma;
        double m_maturity;
        int m_nb_steps_space;
        int m_nb_steps_time;
        double m_theta;
        int m_r;
        double m_dx;
        double m_dt;
        bool vega_computed = false; //boolean pour savoir si on a déjà rempli la grille bumped vol
        coef_eq::CoefEquation* m_alpha;
        coef_eq::CoefEquation* m_beta;
        coef_eq::CoefEquation* m_gamma;
        coef_eq::CoefEquation* m_delta;
        // change to pointers
        payoff::Payoff* m_pf;
        boundary::BoundaryCondition* m_bound_small;
        boundary::BoundaryCondition* m_bound_big;
        std::vector<std::vector<double>> grid_res;
        std::vector<std::vector<double>> grid_res_bumped_sigma;

        double get_price(std::vector<std::vector<double>> grid);
        std::vector<std::vector<double>> main_run();
        std::vector<double> initiate_spot_values(double S0, double sigma, double maturity, int nb_steps_space);


    };

}


#endif // MESH_H_INCLUDED
