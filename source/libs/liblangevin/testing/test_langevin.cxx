#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "../langevin/langevin.hpp"

int main() {

    // conditions initiales et paramètres du système
    langevin::langevin mon_equation; // sac
    mon_equation.x = 3.0;
    mon_equation.theta = 1./5.;
    mon_equation.Q0 = 1.;
    mon_equation.sigma = std::sqrt(2.*mon_equation.x*mon_equation.theta);

    // initialisation
    double Q;
    double V;
    Q = mon_equation.Q0;
    V = mon_equation.define_V0(1.);

    // paramètres des boucles
    //     std::size_t n = 1000;
    double tau; 
    double tau_f;
    double h;
    tau = 0.;
    tau_f =  20.;
    h = 0.001;

    long int seed = 132435;
    srand48(seed);

    std::cout << tau << " " << Q << std::endl;

    //     for (std::size_t i=0 ; i<n ; i++) {
    //          seed += 1;
    //          srand48(seed);

    for (tau = h; tau < tau_f; tau += h) {
        double u1 = drand48();
        double u2 = drand48();
        double rg = std::sqrt(-2.*std::log(u1))*std::cos(2.*M_PI*u2);
        double R = rg * mon_equation.sigma;
        V = mon_equation.eval_V(Q,V,h,R);
        Q = mon_equation.eval_Q(Q,V,h);
        std::cout << tau << " " << Q << std::endl;
    }
    //     }

    return 0;
}
