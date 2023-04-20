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

#include "soil_vangenuchten.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

// Compute soil surface moisture
double VanGenuchten::surfaceWaterContent(const double psi) {
    
    double b        = properties[0]->b;
    double psi_sat  = properties[0]->psi_sat;
    double porosity = properties[0]->porosity;
    double residual = properties[0]->residual;
    double soil_e   = porosity-residual;
    double m        = 1 / (1+b);
    double soil_q   = residual + soil_e*std::pow((1 + std::pow(psi/psi_sat, (1/(1-m)))),-m);
    
    return soil_q;
}

// Estimate soil surface moisture from surface mixing ratio
double VanGenuchten::surfaceWaterContentEstimate(const double sfc_T, const double sfc_q,
                                                 const double atm_p) {
    
    double b        = properties[0]->b;
    double psi_sat  = properties[0]->psi_sat;
    double porosity = properties[0]->porosity;
    double residual = properties[0]->residual;
    double es       = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    double qs       = 0.622*(es/(atm_p-0.378*es));
    double ln       = std::log(sfc_q/qs);
    double soil_e   = porosity-residual;
    double m        = 1 / (1+b);
    double soil_q   = residual+soil_e*std::pow((1+std::pow(c::Rv*sfc_T*ln/(c::grav*psi_sat),1/(1-m))),-m);
    
    return soil_q;
}

// Compute soil water potential (single level)
double VanGenuchten::waterPotential(const double soil_q, const int level) {
    
    double b        = properties[level]->b;
    double psi_sat  = properties[level]->psi_sat;
    double porosity = properties[level]->porosity;
    double residual = properties[level]->residual;
    double Se       = (soil_q-residual)/(porosity-residual);
    double m        = 1 / (1+b);
    double psi      = psi_sat*std::pow((std::pow(Se,-1/m)-1), 1-m);
    
    if (psi>psi_sat) psi = psi_sat;
    
    return psi;
}

// Computes soil moisture conductivity.
double VanGenuchten::conductivityMoisture(const double soil_q, const int level) {

    double b            = properties[level]->b;
    double porosity     = properties[level]->porosity;
    double residual     = properties[level]->residual;
    double K_sat        = properties[level]->K_sat;
    double Se           = (soil_q-residual)/(porosity-residual);
    double m            = 1 / (1+b);
    double conductivity = K_sat*std::sqrt(Se)*std::pow(1 - std::pow(1 - std::pow(Se,1/m),m),2);

    return conductivity;
}

// Computes soil moisture diffusivity
double VanGenuchten::diffusivityMoisture(const double soil_q, const int level) {

    double b            = properties[level]->b;
    double psi_sat      = properties[level]->psi_sat;
    double porosity     = properties[level]->porosity;
    double residual     = properties[level]->residual;
    double K_sat        = properties[level]->K_sat;
    double Se           = (soil_q-residual)/(porosity-residual);
    double soil_e       = porosity-residual;
    double m            = 1 / (1+b);
    double A            = (1-m)*K_sat*psi_sat / (m*soil_e);
    double C            = std::pow(Se,0.5-1/m)*(std::pow(1-std::pow(Se,1/m),-m) +
                          std::pow(1-std::pow(Se,1/m),m) - 2);
    double diffusivity  = A*C;
    
    //  -(1-m)*K_sat*psi_sat*std::pow(Se,-0.5-1/m)*
                        //    std::pow(std::pow(Se,-1/m)-1,-m)*std::pow(1-std::pow(1-std::pow(Se,1/m),m),2) /
                        //    (m*(porosity-residual));

    return diffusivity;
}