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

#include "soil.hpp"
#include "soil_brookscorey.hpp"
#include "soil_campbell.hpp"
#include "soil_vangenuchten.hpp"

#include <cmath>
#include <iostream>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

Soil::Soil(const std::vector<int>& soil_type, const int soil_param, const int soil_model, const int levels) {
    
    // Set up soil properties in the column
    properties.resize(levels);
    for (int i=0; i<levels; i++) {
        properties[i] = SoilType::getProperties(soil_type[i],soil_param);
    }
};

// Compute heat capacity
double Soil::heatCapacity(const double soil_q, const int level) {

    double porosity = properties[level]->porosity;
    double Ci = properties[level]->ci;

    return (1-porosity)*Ci + soil_q*c::Ci_wat + (porosity-soil_q)*c::Cp_air;
}

// Compute surface mixing ratio
double Soil::surfaceMixingRatio(const double sfc_T, const double sfc_q,
                                const double atm_p) {
    
    double psi      = waterPotential(sfc_q, 0);
    double h        = std::exp(c::grav*psi/(c::Rv*sfc_T));
    double es       = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    double hum_sat  = 0.622*(es/(atm_p-0.378*es));
    double hum_spec = h*hum_sat;
    return hum_spec;
}

// Compute soil thermal conductivity
double Soil::conductivityThermal(const double soil_q, const int level) {

    double conductivity;
    double psi = 100.*waterPotential(soil_q,level);
    double pf = std::log10(std::abs(psi));
    if (pf <= 5.1) {
        conductivity = 418.46*std::exp(-(pf+2.7));
    } else {
        conductivity = 0.172;
    }

    return conductivity;
}

// Compute soil thermal diffusivity
double Soil::diffusivityThermal(const double soil_q, const int level) {

    double heat_cap     = heatCapacity(soil_q,level);
    double conductivity = conductivityThermal(soil_q, level);
    double diffusivity  = conductivity / heat_cap;

    return diffusivity;
}

// Factory method that returns the correct soil model
Soil* Soil::getModel(const std::vector<int>& soil_type, const int soil_param, 
                     const int soil_model, const int levels) {
    if (soil_param==1) {
        std::cout<<"[Soil] \t\t Using the Clapp and Hornberger dataset"<<std::endl;
    } else if (soil_param==2) {
        std::cout<<"[Soil] \t\t Using the Cosby et al. dataset"<<std::endl;
    } else if (soil_param==3) {
        std::cout<<"[Soil] \t\t Using the Rawls and Brakensiek dataset"<<std::endl;
    }

    if (soil_model==1) {
        std::cout<<"[Soil] \t\t Using the Brooks-Corey model"<<std::endl;
        return new BrooksCorey(soil_type, soil_param, soil_model, levels);
    } else if (soil_model==2) {
        std::cout<<"[Soil] \t\t Using the Campbell model"<<std::endl;
        return new Campbell(soil_type, soil_param, soil_model, levels);
    } else if (soil_model==3) {
        std::cout<<"[Soil] \t\t Using the van Genuchten model"<<std::endl;
        return new VanGenuchten(soil_type, soil_param, soil_model, levels);
    } else {
        std::cout<<"[Soil] \t\t Invalid soil model: must be an integer 1, 2, or 3"<<std::endl;
        throw(1);
    }
}
