//
//  soil.hpp
//  
//  This namespace handles functions related to soil 
//
//  Created by Jeremy Gibbs on 10/31/17.
//
#include <cmath> 
#include "soil.hpp"
#include "constants.hpp"
#include <iostream>

namespace {
    namespace c = Constants;
}

namespace soil {
    
    // compute surface mixing ratio
    double surfaceMixingRatio(const double psi_nsat, const double porosity, 
                              const double b, const double sfc_T, 
                              const double sfc_q, const double atm_p) {
        
        double psi_n    = -psi_nsat*std::pow((porosity/sfc_q),b);
        double h        = std::exp(c::grav*psi_n/(c::Rv*sfc_T));
        double e        = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
        double hum_sat  = 0.622*(e/(atm_p-0.378*e));
        double hum_spec = h*hum_sat;
        
        //return hum_spec/(1-hum_spec);
        return hum_spec;
    }
    
    // compute soil surface moisture from surface mixing ratio
    double surfaceSoilMoisture(const double psi_nsat, const double porosity, 
                              const double b, const double sfc_T, 
                              const double sfc_r, const double atm_p) {
        
        double es     = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
        double qs     = 0.622*(es/(atm_p-0.378*es));
        double ln     = std::log(sfc_r / (qs*(sfc_r+1)));
        
        return porosity*std::pow((c::grav*psi_nsat)/(c::Rv*sfc_T*ln),(1./b));
    }
    
    // compute average soil thermal conductivity/diffusivity
    std::vector<double> soilThermalTransfer(const std::vector<double> &psi_nsat, 
                                            const std::vector<double> &porosity, 
                                            const std::vector<double> &soil_q, 
                                            const std::vector<double> &b, 
                                            const std::vector<double> &Ci, 
                                            const int depth, const int flag) {
        
        std::vector<double> transfer(depth);
        double psi_n, pf;
        
        // loop through each depth
        for (int d=0; d<depth; ++d) {
            psi_n = 100.*psi_nsat[d]*std::pow((porosity[d]/soil_q[d]),b[d]);
            pf = std::log10(std::abs(psi_n));
            if (pf <= 5.1) {
                transfer[d] = 418.46*std::exp(-(pf+2.7));
            } else {
                transfer[d] = 0.172;
            }
                        
            // convert to thermal diffusivity if flag==1
            if (flag==1) {
                double heat_cap = (1-porosity[d])*Ci[d] + soil_q[d]*c::Ci_wat + (porosity[d]-soil_q[d])*c::Cp_air;
                transfer[d] = transfer[d] / heat_cap;
            }
        }
        
        return transfer;
    }
    
    // compute average soil moisture transfer
    soilTransfer soilMoistureTransfer(const std::vector<double>& psi_nsat, 
                                      const std::vector<double>& K_nsat, 
                                      const std::vector<double>& porosity, 
                                      const std::vector<double>& soil_q, 
                                      const std::vector<double>& b, 
                                      const int depth) {
        
        // struct to hold transfer coefficients
        soilTransfer transfer;
        transfer.d.resize(depth);
        transfer.k.resize(depth);
        
        // loop through each depth
        for (int d=0; d<depth; ++d) {
            transfer.d[d] = -(b[d]*K_nsat[d]*std::abs(psi_nsat[d])/soil_q[d])
                            *std::pow(soil_q[d]/porosity[d],(b[d]+3.));
            transfer.k[d] = K_nsat[d]*std::pow(soil_q[d]/porosity[d],(2.*b[d]+3.));
        }
                
        return transfer;
    }
    
    // set soil type properties at each depth
    // based on Cosby et al. (1984)
    soilProperties soilTypeProperties(const std::vector<int>& soil_type, const int depth) {
        
        // struct to hold soil properties
        soilProperties properties;
        properties.b.resize(depth);
        properties.psi_sat.resize(depth);
        properties.porosity.resize(depth);
        properties.K_sat.resize(depth);
        properties.Ci.resize(depth);
        
        // exponent (unitless)
        const std::vector<double>b_list   = { 2.79,  4.26,  4.74,
                                              5.33,  5.25,  6.77,
                                              8.72,  8.17, 10.73,
                                             10.39, 10.55,  7.75};
        
        // saturation moisture potential (m)
        const std::vector<double>psi_list = {0.02316, 0.01751, 0.03158,
                                             0.06554, 0.04711, 0.03096,
                                             0.05989, 0.04137, 0.02691,
                                             0.04527, 0.05312, 0.35600};
        
        // porosity (volume/volume)
        const std::vector<double>por_list = {0.339, 0.421, 0.434,
                                             0.476, 0.439, 0.404,
                                             0.464, 0.465, 0.406,
                                             0.468, 0.468, 0.863};
        
        // hydraulic conductivity (m/s)
        const std::vector<double>K_list   = {1.60e-5, 9.52e-6, 6.19e-6,
                                             4.73e-6, 5.12e-6, 5.78e-6,
                                             4.11e-6, 4.45e-6, 7.12e-5,
                                             3.43e-6, 2.99e-6, 8.00e-6};
        
        // volumetric heat capacity (J/m^3/K)
        const std::vector<double>Ci_list  = {1.47e6, 1.41e6, 1.34e6,
                                             1.27e6, 1.21e6, 1.18e6,
                                             1.32e6, 1.23e6, 1.18e6,
                                             1.15e6, 1.09e6, 0.84e6};
         
        // loop through each depth to assign soil type properties
        int soil;
        for (int d=0; d<depth; ++d) {
	        
	        soil = soil_type[d] - 1;
	        
            properties.b[d]        = b_list[soil];
            properties.psi_sat[d]  = psi_list[soil];
            properties.porosity[d] = por_list[soil];
            properties.K_sat[d]    = K_list[soil];
            properties.Ci[d]       = Ci_list[soil];
        }
                           
        return properties;
    }  
};
