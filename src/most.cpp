//
//  most.hpp
//  
//  This namespace handles Monin-Obukhov Similarity Theory
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include "most.hpp"
#include "constants.hpp"
#include <cmath>
#include <iostream>

namespace {
    namespace c = Constants;
}

namespace most {
    
    // gradient momentum functions
    double phim(const double zeta) {
        return (zeta >= 0.) ? phimStable(zeta) : phimUnstable(zeta);
    }
    double phimUnstable(const double zeta) {
        return std::pow( (1.-(16.*zeta)), -0.25);
    }
    double phimStable(const double zeta) {
        return 1. + 5.*zeta;
    }
    
    // gradient scalar functions
    double phih(const double zeta) {
        return (zeta >= 0.) ? phihStable(zeta) : phihUnstable(zeta);
    }
    double phihUnstable(const double zeta) {
        return std::pow( (1.-(16.*zeta)), -0.5);
    }
    double phihStable(const double zeta) {
        return 1. + 5.*zeta;
    }

    // integral momentum functions
    double psim(const double zeta) {
        return (zeta >= 0.) ? psimStable(zeta) : psimUnstable(zeta);
    }
    double psimStable(const double zeta) {
        return -5.*zeta;
    }
    double psimUnstable(const double zeta) {
        double x = std::pow( (1.-(16.*zeta)), 0.25);
        return 2.*std::log((1.+x)/2.)+std::log((1+std::pow(x,2.))/2.)-2*std::atan(x)+c::pi/2.;
    }
    
    // integral scalar functions
    double psih(const double zeta) {
        return (zeta >= 0.) ? psihStable(zeta) : psihUnstable(zeta);
    }
    double psihStable(const double zeta) {
        return -5.*zeta;
    }
    double psihUnstable(const double zeta) {
        double x = std::pow( (1.-(16.*zeta)), 0.25);
        return 2.*std::log(1+std::pow(x,2)/2.);
    }
    
    // log-law functions
    double fm(const double zm_over_zo, const double zeta_m, const double zeta_o) {
        return c::vonk / (std::log(zm_over_zo) - psim(zeta_m) + psim(zeta_o));
    }
    double fh(const double zs_over_zt, const double zeta_s, const double zeta_t) {
        return c::vonk / (std::log(zs_over_zt) - psih(zeta_s) + psih(zeta_t));
    }
};
