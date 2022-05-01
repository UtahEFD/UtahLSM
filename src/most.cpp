/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2022 Jeremy A. Gibbs
 * Copyright (c) 2017–2022 Rob Stoll
 * Copyright (c) 2017–2022 Eric Pardyjak
 * Copyright (c) 2017–2022 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "most.hpp"

#include <cmath>
#include <iostream>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

namespace most {
    
    // Gradient momentum function
    double phim(const double zeta) {
        return (zeta >= 0.) ? phimStable(zeta) : phimUnstable(zeta);
    }

    // Gradient momentum function under neutral/stable conditions
    double phimStable(const double zeta) {
        return 1. + 5.*zeta;
    }

    // Gradient momentum function under unstable conditions
    double phimUnstable(const double zeta) {
        return std::pow( (1.-(16.*zeta)), -0.25);
    }
    
    // Gradient momentum function
    double phih(const double zeta) {
        return (zeta >= 0.) ? phihStable(zeta) : phihUnstable(zeta);
    }
    // Gradient scalar function under neutral/stable conditions
    double phihStable(const double zeta) {
        return 1. + 5.*zeta;
    }

    // Gradient scalar function under unstable conditions
    double phihUnstable(const double zeta) {
        return std::pow( (1.-(16.*zeta)), -0.5);
    }

    // Integral stability correction for momentum
    double psim(const double zeta) {
        return (zeta >= 0.) ? psimStable(zeta) : psimUnstable(zeta);
    }

    // Integral stability correction for momentum under neutral/stable conditions
    double psimStable(const double zeta) {
        return -5.*zeta;
    }

    // Integral stability correction for momentum under unstable conditions
    double psimUnstable(const double zeta) {
        double x = std::pow( (1.-(16.*zeta)), 0.25);
        return 2.*std::log((1.+x)/2.)+std::log((1+std::pow(x,2.))/2.)-2*std::atan(x)+c::pi/2.;
    }
    
    // Integral stability correction for scalars
    double psih(const double zeta) {
        return (zeta >= 0.) ? psihStable(zeta) : psihUnstable(zeta);
    }

    // Integral stability correction for scalars under neutral/stable conditions
    double psihStable(const double zeta) {
        return -5.*zeta;
    }

    // Integral stability correction for scalars under unstable conditions
    double psihUnstable(const double zeta) {
        double y = std::pow( (1.-(16.*zeta)), 0.5);
        return 2.*std::log((1+y)/2.);
    }
    
    // Common log-law function for momentum
    double fm(const double zm_over_zo, const double zeta_m, const double zeta_o) {
        return c::vonk / (std::log(zm_over_zo) - psim(zeta_m) + psim(zeta_o));
    }

    // Common log-law function for heat
    double fh(const double zs_over_zt, const double zeta_s, const double zeta_t) {
        return c::vonk / (std::log(zs_over_zt) - psih(zeta_s) + psih(zeta_t));
    }
};