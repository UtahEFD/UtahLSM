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
    soilProperties soilTypeProperties(const std::vector<int>& soil_type, const int depth) {
        
        // struct to hold soil properties
        soilProperties properties;
        properties.b.resize(depth);
        properties.psi_sat.resize(depth);
        properties.porosity.resize(depth);
        properties.K_sat.resize(depth);
        properties.Ci.resize(depth);
        
        
        
        // exponent (unitless)
        const std::vector<double>b_list   = {4.05, 4.38,  4.90,  5.30,  5.39, 7.12, 
                                             7.75, 8.52, 10.40, 10.40, 11.40, 7.75};
        // saturation moisture potential (cm)
        //const std::vector<double>psi_list = {-12.1,  -9.0, -21.8, -78.6, -47.8, -29.9,
        //                                     -35.6, -63.0, -15.3, -49.0, -40.5, -35.6};
        
        const std::vector<double>psi_list = {12.1,  9.0, 21.8, 78.6, 47.8, 29.9,
                                             35.6, 63.0, 15.3, 49.0, 40.5, 35.6};
        
        // porosity (unitless)
        const std::vector<double>por_list = {0.395, 0.410, 0.435, 0.485, 0.451, 0.420,
                                             0.477, 0.476, 0.426, 0.492, 0.482, 0.863};
        // hydraulic conductivity (cm/s)
        const std::vector<double>K_list   = {.01760, .01563, .00341, .00072, .00070, .00063,
                                             .00017, .00025, .00022, .00010, .00013, .00080};
        // volumetric heat capacity (J/cm^3/K)
        const std::vector<double>Ci_list  = {1.47, 1.41, 1.34, 1.27, 1.21, 1.18, 
                                             1.32, 1.23, 1.18, 1.15, 1.09, 0.84};
         
        // loop through each depth to assign soil type properties
        int soil;
        for (int d=0; d<depth; ++d) {
	        
	        soil = soil_type[d] - 1;
	        
            properties.b[d]        = b_list[soil];             // unitless
            properties.psi_sat[d]  = psi_list[soil] / 100.;    // from cm to m
            properties.porosity[d] = por_list[soil];           // unitless
            properties.K_sat[d]    = K_list[soil] / 100.;      // from cm/s to m/s
            properties.Ci[d]       = Ci_list[soil] * 1000000.; // from J/cm^3/K to J/m^3/K
        }
                           
        return properties;
    }  
};
