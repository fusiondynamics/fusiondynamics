#include <cmath>

#include "langevin.hpp"

namespace langevin { 

    double langevin::define_V0(const double f_) const {
        double V0;
        V0 = -(x + std::sqrt(1. + x * x));
        return V0;
    }

    double langevin::eval_V(const double Q_, const double V_, const double h_, const double R_) const{
        double V;
        V = V_*(1. - 2.*x*h_) + Q_*h_ + R_ * h_;
        return V;
    }

    double langevin::eval_Q(const double Q_, const double V_, const double h_) const {
        double Q;
        Q = Q_ + V_*h_;
        return Q;
    }

} // namespace langevin
