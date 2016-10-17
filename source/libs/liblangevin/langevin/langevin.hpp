#ifndef LANGEVIN__LANGEVIN_HPP
#define LANGEVIN__LANGEVIN_HPP

namespace langevin {

    struct langevin {

        // methodes
        double define_V0(const double f_) const;
        double eval_V(const double Q_, const double V_, const double h_, const double R_) const;
        double eval_Q(const double Q_, const double V_, const double h_) const;

        // attributs
        double x;
        double theta;
        double sigma;
        double Q0;

    };

} // namespace langevin

#endif // langevin_HPP
