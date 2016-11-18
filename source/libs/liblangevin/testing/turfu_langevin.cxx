#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>

#include "../langevin/langevin.hpp"


int main() {


    //***********************************************//
    // conditions initiales et paramètres du système //
    //***********************************************//

    langevin::langevin mon_equation; // sac
    mon_equation.x = 3.0;

    // on demande à l'utilisateur la valeur de theta
    std::clog << "theta = ";
    std::cin >> mon_equation.theta;

    mon_equation.Q0 = 1.;
    mon_equation.sigma = std::sqrt(2.*mon_equation.x*mon_equation.theta);


    //****************//
    // initialisation //
    //****************//

    double Q;
    double V;

    // on demande à l'utilisateur la valeur de K/Beff
    double f;
    std::clog << "K/Beff = ";
    std::cin >> f;


    //*********//
    // routine //
    //*********//

    // choix du nombre de trajectoires
    long int seed;
    long int last_seed;
    last_seed = 1000;

    // paramètres de la boucle sur le temps
    double tau;
    double tau_f;
    double h;
    tau_f = 20.;
    h = 0.0001;

    std::vector<double> Nb;
    int i = 0;
    for (tau = 0.; tau < tau_f; tau+=h) {
        Nb.push_back(0.);
        i += 1;
    }

    for (seed = 0; seed < last_seed; seed++) {

        i = 0;
        std::clog << 100.*seed/last_seed << "%" << std::endl;
        Q = mon_equation.Q0;
        V = mon_equation.define_V0(f);
        srand48(seed);

        for (tau = h; tau < tau_f; tau += h) {
            i += 1;
            double old_Q = Q; // sauvegarde
            double u1 = drand48();
            double u2 = drand48();
            double rg = std::sqrt(-2.*std::log(u1))*std::cos(2.*M_PI*u2);
            double R = rg * mon_equation.sigma;
            V = mon_equation.eval_V(Q,V,h,R);
            Q = mon_equation.eval_Q(Q,V,h);

            if (Q < 0. && old_Q > 0.) {Nb[i] = Nb[i] + 1.;}
            if (Q > 0. && old_Q < 0.) {Nb[i] = Nb[i] - 1.;}
        }

//         std::clog << Nb[150000] << std::endl;

    }

    i = 0;
    for(tau = 0; tau < tau_f; tau += h) {
        std::cout << tau << " " << Nb[i]/last_seed << std::endl;
        i += 1;
    }

    return 0;
}
