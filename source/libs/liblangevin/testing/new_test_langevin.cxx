#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>

#include "../langevin/langevin.hpp"

int main() {

    // conditions initiales et paramètres du système
    langevin::langevin mon_equation; // sac
    mon_equation.x = 3.0;
    mon_equation.theta = 1./5.;
    mon_equation.Q0 = 1.;
    mon_equation.sigma = std::sqrt(2.*mon_equation.x*mon_equation.theta);

    // initialisation
    std::vector<double> Q;
    std::vector<double> V;
    long int seed;
    for (seed = 0; seed < 10000; seed++) {
        Q.push_back(mon_equation.Q0);
        V.push_back(mon_equation.define_V0(0.));
    }

    // paramètres de la boucle sur le temps
    double tau;
    double tau_f;
    double h;
    tau = 0.;
    tau_f =  20.;
    h = 0.0001;

    std::cout << tau << " " << 0. << std::endl;

    // début de la routine
    for (tau = h; tau < tau_f; tau += h) {
        double count = 0.;
        double sum_Q = 0.;
        for(seed = 0; seed < 10000; seed++) {
            srand48(seed);
            double u1 = drand48();
            double u2 = drand48();
            double rg = std::sqrt(-2.*std::log(u1))*std::cos(2.*M_PI*u2);
            double R = rg * mon_equation.sigma;
            V[seed] = mon_equation.eval_V(Q[seed],V[seed],h,R);
            Q[seed] = mon_equation.eval_Q(Q[seed],V[seed],h);
            sum_Q += Q[seed];
            if(Q[seed] < 0.) {
                count += 1.;
            }
        }
        std::cout << tau << " " << sum_Q/10000. << " " << count/10000. << std::endl;
    }

    return 0;
}
