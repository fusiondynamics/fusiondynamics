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

    /* création de deux tableaux ou seront stockées les
       valeurs de Q et de V */
    std::vector<double> Q;
    std::vector<double> V;

    // on demande à l'utilisateur la valeur de K/Beff
    double f;
    std::clog << "K/Beff = ";
    std::cin >> f;

    // choix du nombre de trajectoires
    long int seed;
    long int last_seed;
    last_seed = 10000;

    /* initialisation des tableaux avec les conditions
       initiales */
    for (seed = 0; seed < last_seed; seed++) {
        Q.push_back(mon_equation.Q0);
        V.push_back(mon_equation.define_V0(f));
    }


    //*********//
    // routine //
    //*********//

    // paramètres de la boucle sur le temps
    double tau;
    double tau_f;
    double h;
    tau = 0.;
    tau_f =  20.;
    h = 0.0001;

    std::cout << tau << " " << mon_equation.Q0 << " " << 0 << std::endl;

    // début de la routine
    for (tau = h; tau < tau_f; tau += h) { // boucle sur le temps

        // remise à zéro des compteurs
        double count = 0.;
        double sum_Q = 0.;

        for(seed = 0; seed < last_seed; seed++) { // boucle sur les trajectoires

            // méthode de Box-Muller
            srand48(seed);
            double u1 = drand48();
            double u2 = drand48();
            double rg = std::sqrt(-2.*std::log(u1))*std::cos(2.*M_PI*u2);
            double R = rg * mon_equation.sigma;

            // calcul de V puis de Q
            V[seed] = mon_equation.eval_V(Q[seed],V[seed],h,R);
            Q[seed] = mon_equation.eval_Q(Q[seed],V[seed],h);

            // on somme la valeur de Q à l'instant tau de chacune des trajectoires
            sum_Q += Q[seed];

            /* si Q est négatif pour la trajectoire considéré à l'instaut, on ajoute
               1 au compteur */
            if(Q[seed] < 0.) {
                count += 1.;
            }

        }

        // affichage sur 3 colonnes : tau, Q(tau), P(tau)
        std::cout << tau << " " << sum_Q/last_seed << " " << count/last_seed << std::endl;

    }

    return 0;
}
