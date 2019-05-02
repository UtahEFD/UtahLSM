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

#include <cmath>
#include <iostream>

#include "constants.hpp"

namespace {
    namespace c = constants;
}

Soil::Soil(const std::vector<int>& soil_type, const int soil_param, const int soil_model, const int levels) {
               
    // Set up soil properties in the column
    for (int i=0; i<levels; i++) {
        properties[i] = SoilType::get(soil_type[i],soil_param);
    }
};

// Compute heat capacity
double Soil::heatCapacity(const double porosity, const double Ci, const double soil_q) {
    return (1-porosity)*Ci + soil_q*c::Ci_wat + (porosity-soil_q)*c::Cp_air;
}

// Compute surface mixing ratio
double Soil::surfaceMixingRatio(const double psi_sat, const double porosity,
                          const double residual,  const double b,
                          const double sfc_T, const double sfc_q,
                          const double atm_p, const int model) {
    double psi      = waterPotential(psi_sat, porosity, residual, sfc_q, b, model);
    double h        = std::exp(c::grav*psi/(c::Rv*sfc_T));
    double es        = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    double hum_sat  = 0.622*(es/(atm_p-0.378*es));
    double hum_spec = h*hum_sat;
    return hum_spec;
}

// Compute soil surface moisture
double Soil::surfaceWaterContent(const double psi, const double psi_sat,
                           const double porosity, const double residual,
                           const double b, const int model) {
    
    double soil_q;
    double soil_e = porosity-residual;
    
    // model=1: Brooks and Corey (1964)
    // model=2: Campbell (1974)
    // model=3: van Genuchten (1980)
    if (model==1) {
        soil_q = residual+soil_e*std::pow(psi_sat/psi,(1./b));
    } else if (model==2) {
        soil_q = porosity*std::pow(psi_sat/psi,(1./b));
    } else if (model==3) {
        double m = 1 / (1+b);
        soil_q = residual+soil_e*std::pow( (1 / (1 + std::pow(psi/psi_sat, (1/(1-m))))),m);
    } else {
        std::cout<<"Soil model must be 1, 2, or 3"<<std::endl;
        throw(1);
    }
    return soil_q;
    
}

// Estimate soil surface moisture from surface mixing ratio
double Soil::surfaceWaterContentEstimate(const double psi_sat, const double porosity,
                           const double residual, const double b,
                           const double sfc_T, const double sfc_q,
                           const double atm_p, const int model) {
    
    double soil_q;
    double es     = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
    double qs     = 0.622*(es/(atm_p-0.378*es));
    double ln     = std::log(sfc_q/qs);
    double soil_e = porosity-residual;
    
    // model=1: Brooks and Corey (1964)
    // model=2: Campbell (1974)
    // model=3: van Genuchten (1980)
    if (model==1) {
        soil_q = residual+soil_e*std::pow(c::Rv*sfc_T*ln/(c::grav*psi_sat),-1./b);
    } else if (model==2) {
        soil_q = porosity*std::pow(c::Rv*sfc_T*ln/(c::grav*psi_sat),-1./b);
    } else if (model==3) {
        double m = 1 / (1+b);
        soil_q = residual+soil_e*std::pow(1/ (1+std::pow(c::Rv*sfc_T*ln/(c::grav*psi_sat),1/(1-m))),m);
    } else {
        std::cout<<"Soil model must be 1, 2, or 3"<<std::endl;
        throw(1);
    }
    return soil_q;
}

// Compute soil thermal conductivity/diffusivity
struct Soil::ThermalTransfer Soil::thermalTransfer(const std::vector<double> &psi_sat,
                                                   const std::vector<double> &porosity,
                                                   const std::vector<double> &residual,
                                                   const std::vector<double> &soil_q,
                                                   const std::vector<double> &b,
                                                   const std::vector<double> &Ci,
                                                   const int depth, const int model) {
    
    struct ThermalTransfer transfer;
    transfer.d.resize(depth);
    transfer.k.resize(depth);
    double psi, pf;
    
    // loop through each depth
    for (int d=0; d<depth; ++d) {
        
        psi = 100.*waterPotential(psi_sat[d],porosity[d],residual[d],soil_q[d],b[d],model);
        //)  psi_nsat[d]*std::pow((porosity[d]/soil_q[d]),b[d]);
        pf = std::log10(std::abs(psi));
        if (pf <= 5.1) {
            transfer.k[d] = 418.46*std::exp(-(pf+2.7));
        } else {
            transfer.k[d] = 0.172;
        }
                    
        // convert to thermal diffusivity
        double heat_cap = heatCapacity(porosity[d], Ci[d], soil_q[d]);
        transfer.d[d] = transfer.k[d] / heat_cap;
    }
    
    return transfer;
}

// Compute average soil moisture transfer
struct Soil::MoistureTransfer Soil::moistureTransfer(const std::vector<double>& psi_sat,
                                                     const std::vector<double>& K_sat,
                                                     const std::vector<double>& porosity,
                                                     const std::vector<double>& residual,
                                                     const std::vector<double>& soil_q,
                                                     const std::vector<double>& b,
                                                     const int depth, const int model) {
    
    // Local variables
    double Se;
    struct MoistureTransfer transfer;
    transfer.k.resize(depth);
    transfer.d.resize(depth);
    // Loop through each depth
    for (int d=0; d<depth; ++d) {
        
        Se = (soil_q[d]-residual[d])/(porosity[d]-residual[d]);
        
        // model=1: Brooks and Corey (1964)
        // model=2: Campbell (1974)
        // model=3: van Genuchten (1980)
        if (model==1) {
            transfer.k[d] = K_sat[d]*std::pow(Se,(2.*b[d]+3.));
            transfer.d[d] = -b[d]*K_sat[d]*psi_sat[d]*std::pow(Se,(b[d]+2.))/(porosity[d]-residual[d]);
        } else if (model==2) {
            transfer.k[d] = K_sat[d]*std::pow(soil_q[d]/porosity[d],(2.*b[d]+3.));
            transfer.d[d] = -b[d]*K_sat[d]*psi_sat[d]*std::pow(soil_q[d]/porosity[d],(b[d]+2.))/soil_q[d];
        } else if (model==3) {
            double m = 1 / (1+b[d]);
            transfer.k[d] = K_sat[d]*std::sqrt(Se)*std::pow(1 - std::pow(1 - std::pow(Se,1/m),m),2);
            transfer.d[d] = (1-m)*K_sat[d]*psi_sat[d]*std::pow(Se,0.5-1/m)*
            (std::pow(1-std::pow(Se,1/m),-m) + std::pow(1-std::pow(Se,1/m),m) -2)/(m*(porosity[d]-residual[d]));
        }
    }
            
    return transfer;
}

// Compute soil water potential (single level)
double Soil::waterPotential(const double psi_sat, const double porosity,
                      const double residual, const double soil_q,
                      const double b, const int model) {
    
    // Local variables
    double psi, Se;
    
    Se = (soil_q-residual)/(porosity-residual);
    
    // model=1: Brooks and Corey (1964)
    // model=2: Campbell (1974)
    // model=3: van Genuchten (1980)
    if (model==1) {
        psi = psi_sat*std::pow(Se,-b);
        if (psi>psi_sat) psi = psi_sat;
    } else if (model==2) {
        psi = psi_sat*std::pow(soil_q/porosity,-b);
        if (psi>psi_sat) psi = psi_sat;
    } else if (model==3) {
        double m = 1 / (1+b);
        psi = psi_sat*std::pow((std::pow(Se,-1/m)-1), 1-m);
        if (psi>psi_sat) psi = psi_sat;
    } else {
        std::cout<<"Soil model must be 1, 2, or 3"<<std::endl;
        throw(1);
    }
    return psi;
}

// Compute soil water potential (full column)
std::vector<double> Soil::waterPotential(const std::vector<double>& psi_sat,
                                   const std::vector<double>& porosity,
                                   const std::vector<double>& residual,
                                   const std::vector<double>& soil_q,
                                   const std::vector<double>& b,
                                   const int depth, const int model) {
    
    // Local variables
    double Se;
    std::vector<double> psi(depth);
    
    // Loop through each depth
    for (int d=0; d<depth; ++d) {
        
        Se = (soil_q[d]-residual[d])/(porosity[d]-residual[d]);
        
        // model=1: Brooks and Corey (1964)
        // model=2: Campbell (1974)
        // model=3: van Genuchten (1980)
        if (model==1) {
            psi[d] = psi_sat[d]*std::pow(Se,-b[d]);
            if (psi[d]>psi_sat[d]) psi[d] = psi_sat[d];
        } else if (model==2) {
            psi[d] = psi_sat[d]*std::pow(soil_q[d]/porosity[d],-b[d]);
            if (psi[d]>psi_sat[d]) psi[d] = psi_sat[d];
        } else if (model==3) {
            double m = 1 / (1+b[d]);
            psi[d] = psi_sat[d]*std::pow((std::pow(Se,-1/m)-1), 1-m);
            if (psi[d]>psi_sat[d]) psi[d] = psi_sat[d];
        } else {
            std::cout<<"Soil model must be 1, 2, or 3"<<std::endl;
            throw(1);
        }
    }
    return psi;
}