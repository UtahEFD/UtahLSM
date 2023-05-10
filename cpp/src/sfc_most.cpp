/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2023 Jeremy A. Gibbs
 * Copyright (c) 2017–2023 Rob Stoll
 * Copyright (c) 2017–2023 Eric Pardyjak
 * Copyright (c) 2017–2023 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "sfc_most.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

// Gradient momentum function
double MOST::phim(const double z, const double obukL) {
    double zeta = (obukL == 0) ? 0 : z/obukL; 
    return (zeta >= 0.) ? phimStable(zeta) : phimUnstable(zeta);
}

// Gradient momentum function under neutral/stable conditions
double MOST::phimStable(const double zeta) {
    return 1. + 5.*zeta;
}

// Gradient momentum function under unstable conditions
double MOST::phimUnstable(const double zeta) {
    return std::pow( (1.-(16.*zeta)), -0.25);
}

// Gradient momentum function
double MOST::phih(const double z, const double obukL) {
    double zeta = (obukL == 0) ? 0 : z/obukL; 
    return (zeta >= 0.) ? phihStable(zeta) : phihUnstable(zeta);
}

// Gradient scalar function under neutral/stable conditions
double MOST::phihStable(const double zeta) {
    return 1. + 5.*zeta;
}

// Gradient scalar function under unstable conditions
double MOST::phihUnstable(const double zeta) {
    return std::pow( (1.-(16.*zeta)), -0.5);
}

// Integral stability correction for momentum
double MOST::psim(const double z, const double obukL) {
    double zeta = (obukL == 0) ? 0 : z/obukL;
    return (zeta >= 0.) ? psimStable(zeta) : psimUnstable(zeta);
}

// Integral stability correction for momentum under neutral/stable conditions
double MOST::psimStable(const double zeta) {
    return -5.*zeta;
}

// Integral stability correction for momentum under unstable conditions
double MOST::psimUnstable(const double zeta) {
    
    double x = std::pow( (1.-(16.*zeta)), 0.25);
    double log1 = 2.*std::log((1.+x)/2.);
    double log2 = std::log((1.+std::pow(x,2.))/2.);
    double atan = 2.* std::cos(x) / std::sin(x);
    
    //double atan 2.*std::atan(x);
    double pio2 = c::pi/2.;
    
    std::cout<<std::defaultfloat;
    std::cout<<std::setprecision(17);
    std::cout<<"PSIMU\t\t"<<"zeta"<<": "<<zeta<<std::endl;
    std::cout<<"PSIMU\t\t"<<"x"<<": "<<x<<std::endl;
    std::cout<<"PSIMU\t\t"<<"log1"<<": "<<log1<<std::endl;
    std::cout<<"PSIMU\t\t"<<"log2"<<": "<<log2<<std::endl;
    std::cout<<"PSIMU\t\t"<<"atan"<<": "<<atan<<std::endl;
    std::cout<<"PSIMU\t\t"<<"pio2"<<": "<<pio2<<std::endl;
    
    return 2.*std::log((1.+x)/2.)+std::log((1.+std::pow(x,2.))/2.)-2.*std::atan2(1.,std::pow( (1.-(16.*zeta)), -0.25))+c::pi/2.;
}

// Integral stability correction for scalars
double MOST::psih(const double z, const double obukL) {
    double zeta = (obukL == 0) ? 0 : z/obukL;
    return (zeta >= 0.) ? psihStable(zeta) : psihUnstable(zeta);
}

// Integral stability correction for scalars under neutral/stable conditions
double MOST::psihStable(const double zeta) {
    return -5.*zeta;
}

// Integral stability correction for scalars under unstable conditions
double MOST::psihUnstable(const double zeta) {
    double y = std::pow( (1.-(16.*zeta)), 0.5);
    double log1 = 2.*std::log((1.+y)/2.);
    std::cout<<std::defaultfloat;
    std::cout<<std::setprecision(17);
    std::cout<<"PSIHU\t\t"<<"zeta"<<": "<<zeta<<std::endl;
    std::cout<<"PSIHU\t\t"<<"x"<<": "<<y<<std::endl;
    std::cout<<"PSIHU\t\t"<<"log1"<<": "<<log1<<std::endl;
    return 2.*std::log((1.+y)/2.);
}

// Common log-law function for momentum
double MOST::fm(const double z1, const double z0, const double obukL) {
    return c::vonk / (std::log(z1/z0) - psim(z1,obukL) + psim(z0,obukL));
}

// Common log-law function for heat
double MOST::fh(const double z1, const double z0, const double obukL) {
    return c::vonk / (std::log(z1/z0) - psih(z1,obukL) + psih(z0,obukL));
}