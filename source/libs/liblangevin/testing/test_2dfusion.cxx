#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

int main(void) {

    // initialization

    double Q;
    double V;

    double h = 1e-4;
    double x = 3.;
    double theta = 0.5;
    double sigma = std::sqrt(2.*x*theta*h);

    double Q0 = 1.;
    double f = 2.;
    double V0 = -(x + std::sqrt(1. + x * x))*sqrt(f);

    int trajectory;
    int last_trajectory = 10000;
    long int seed = 0;

    double tau;
    double tau_i = 0.;
    double tau_f = 0.4;

    double vector_size = tau_f/h + 1;
    std::vector<double> v_mean_Q(vector_size, 0.0);
    std::vector<double> proba(vector_size, 0.0);
    std::vector<double> v_Q(vector_size, 0.0);
    std::vector<double> v_tau(vector_size, 0.0);


    for (trajectory = 0; trajectory < last_trajectory; trajectory++) {

        // initial conditions
        Q = Q0;
        V = V0;
        v_Q.at(0) = Q;
        tau = tau_i;

        srand48(seed); // seed updating

        for (int j = 0; j < v_tau.size(); j++) {

            v_tau.at(j) = tau;
            tau += h;

            // Box-Muller
            double u1 = drand48();
            double u2 = drand48();
            double rg = std::sqrt(-2.*std::log(u1))*std::cos(2.*M_PI*u2);
            double R = rg * sigma;

            // Euler method
            V = V*(1. - 2.*x*h) + Q*h + R;
            Q = Q + V*h;
            v_Q.at(j) = Q;

        } // time loop

        std::clog << trajectory << std::endl;

        // mean trajectory and fusion probability computation loop
        for(int k = 0; k < proba.size(); k++) {
            v_mean_Q.at(k) = v_mean_Q.at(k) + v_Q.at(k);
            if(v_Q.at(k) < 0.) {
                proba.at(k) = proba.at(k) + 1;
            }
        }

        seed += 1; // changing the seed

    } // trajectory loop

    // print results
    for(int l = 0; l < v_mean_Q.size(); l++) {
        std::cout << v_tau.at(l) << " " << v_mean_Q.at(l)/last_trajectory << " " << proba.at(l)/last_trajectory << std::endl;

    }

    return 0;
}
