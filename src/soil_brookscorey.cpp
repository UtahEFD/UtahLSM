/*
 * UtahLSM
 * 
 * Copyright (c) 2019 Jeremy A. Gibbs
 * Copyright (c) 2019 Pete Willemsen
 * Copyright (c) 2019 Rob Stoll
 * Copyright (c) 2019 Eric Pardyjak
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "soil_brookscorey.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

// Compute soil surface moisture
double BrooksCorey::surfaceWaterContent(const double psi) {
    
    double b        = properties[0]->b;
    double psi_sat  = properties[0]->psi_sat;
    double porosity = properties[0]->porosity;
    double residual = properties[0]->residual;
    double soil_e   = porosity-residual;
    double soil_q   = residual+soil_e*std::pow(psi_sat/psi,1./b);
    //std::cout<<"SFW: "<<std::setprecision(8)<<"soilq: "<<soil_q<<" resid: "<<residual<<" porosity: "<<porosity<<" psi: "<<psi<<" psisat: "<<psi_sat<<std::endl;
    return soil_q;
}

// Estimate soil surface moisture from surface mixing ratio
double BrooksCorey::surfaceWaterContentEstimate(const double sfc_T, const double sfc_q,
                                                const double atm_p) {
    
    double b        = properties[0]->b;
    double psi_sat  = properties[0]->psi_sat;
    double porosity = properties[0]->porosity;
    double residual = properties[0]->residual;
    double es       = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    double qs       = 0.622*(es/(atm_p-0.378*es));
    double ln       = std::log(sfc_q/qs);
    double soil_e   = porosity-residual;
    double soil_q   = residual+soil_e*std::pow(c::Rv*sfc_T*ln/(c::grav*psi_sat),-1./b);
    //std::cout<<"sfc_q: "<<sfc_q<<" soil_q: "<<soil_q<<std::endl;
    return soil_q;
}

// Compute soil water potential (single level)
double BrooksCorey::waterPotential(const double soil_q, const int level) {
    
    double b        = properties[level]->b;
    double psi_sat  = properties[level]->psi_sat;
    double porosity = properties[level]->porosity;
    double residual = properties[level]->residual;
    double Se       = (soil_q-residual)/(porosity-residual);
    double psi      = psi_sat*std::pow(Se,-b);
    
    double h        = std::exp(c::grav*psi/(c::Rv*295));
    std::cout<<"--- PSI (level "<<level<<"): "<<std::setprecision(4)<<"soilq: "<<soil_q<<" resid: "<<residual<<" porosity: "<<porosity<<" Se: "<<Se<<" psi: "<<psi<<" psisat: "<<psi_sat<<" h: "<<h<<std::endl;
    
    if (psi>psi_sat) {
        std::cout<<"    ====== WARNING: PSI (level "<<level<<"): Saturated! ======"<<std::endl;
        double h        = std::exp(c::grav*psi/(c::Rv*295));
        psi = psi_sat;
        Se = 1;
        std::cout<<"--- FIXED (level "<<level<<"): "<<std::setprecision(4)<<"soilq: "<<soil_q<<" resid: "<<residual<<" porosity: "<<porosity<<" Se: "<<Se<<" psi: "<<psi<<" psisat: "<<psi_sat<<" h: "<<h<<std::endl;
        std::cout<<"    ================================================"<<std::endl;
    }
    
    if (soil_q<0) {
        std::exit(0);
    }
    if (soil_q > 1) {
        std::exit(0);
    }
    return psi;
}

// Computes soil moisture conductivity
double BrooksCorey::conductivityMoisture(const double soil_q, const int level) {

    double b            = properties[level]->b;
    double porosity     = properties[level]->porosity;
    double residual     = properties[level]->residual;
    double K_sat        = properties[level]->K_sat;
    double Se           = (soil_q-residual)/(porosity-residual);
    double conductivity = K_sat*std::pow(Se,(2.*b+3.));
    //std::cout<<"CON (level "<<level<<"): "<<std::setprecision(17)<<"soilq: "<<soil_q<<" resid: "<<residual<<" porosity: "<<porosity<<" Se: "<<Se<<" b "<<b<<" Ksat "<<K_sat<<" K: "<<conductivity<<std::endl;
    return conductivity;
}

// Computes soil moisture diffusivity
double BrooksCorey::diffusivityMoisture(const double soil_q, const int level) {

    double b            = properties[level]->b;
    double psi_sat      = properties[level]->psi_sat;
    double porosity     = properties[level]->porosity;
    double residual     = properties[level]->residual;
    double K_sat        = properties[level]->K_sat;
    double Se           = (soil_q-residual)/(porosity-residual);
    double diffusivity  = -b*K_sat*psi_sat*std::pow(Se,(b+2.))/(porosity-residual);
    //std::cout<<"DIF (level "<<level<<"): "<<std::setprecision(17)<<"soilq: "<<soil_q<<" resid: "<<residual<<" porosity: "<<porosity<<" Se: "<<Se<<" b "<<b<<" D: "<<diffusivity<<std::endl;
    return diffusivity;
}