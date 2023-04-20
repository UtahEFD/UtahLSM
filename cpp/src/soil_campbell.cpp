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

#include "soil_campbell.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

// Compute soil surface moisture
double Campbell::surfaceWaterContent(const double psi) {
    
    double b        = properties[0]->b;
    double psi_sat  = properties[0]->psi_sat;
    double porosity = properties[0]->porosity;
    double soil_q   = porosity*std::pow(std::abs(psi_sat/psi),(1./b));
    
    return soil_q;
}

// Estimate soil surface moisture from surface mixing ratio
double Campbell::surfaceWaterContentEstimate(const double sfc_T, const double sfc_q,
                                             const double atm_p) {
    double b        = properties[0]->b;
    double psi_sat  = properties[0]->psi_sat;
    double porosity = properties[0]->porosity;
    double es       = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    double qs       = 0.622*(es/(atm_p-0.378*es));
    double ln       = std::log(sfc_q/qs);
    double soil_q   = porosity*std::pow(c::Rv*sfc_T*ln/(c::grav*psi_sat),-1./b);
    
    return soil_q;
}

// Compute soil water potential (single level)
double Campbell::waterPotential(const double soil_q, const int level) {
    
    double b        = properties[level]->b;
    double psi_sat  = properties[level]->psi_sat;
    double porosity = properties[level]->porosity;
    double psi      = psi_sat*std::pow(soil_q/porosity,-b);
    if (psi>psi_sat) psi = psi_sat;
    
    return psi;
}

// Computes soil moisture conductivity.
double Campbell::conductivityMoisture(const double soil_q, const int level) {

    double b            = properties[level]->b;
    double porosity     = properties[level]->porosity;
    double K_sat        = properties[level]->K_sat;
    double conductivity = K_sat*std::pow(soil_q/porosity,(2.*b+3.));

    return conductivity;
}

// Computes soil moisture diffusivity
double Campbell::diffusivityMoisture(const double soil_q, const int level) {

    double b            = properties[level]->b;
    double psi_sat      = properties[level]->psi_sat;
    double porosity     = properties[level]->porosity;
    double K_sat        = properties[level]->K_sat;
    double diffusivity  = -b*K_sat*psi_sat*std::pow(soil_q/porosity,(b+2.))/porosity;

    return diffusivity;
}