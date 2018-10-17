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
        
        double psi_n    = psi_nsat*std::pow((porosity/sfc_q),b);
        double h        = std::exp(c::grav*psi_n/(c::Rv*sfc_T));
        //double w        = sfc_q/porosity;
        //double h        = w < 0.75 ? w/0.75 : 1;
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
        double ln     = std::log(sfc_r / qs);
        
        return porosity*std::pow((std::abs(c::grav*psi_nsat))/(c::Rv*sfc_T*std::abs(ln)),(1./b));
    }
    
    // compute soil thermal conductivity/diffusivity
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
        
        // struct to hold transfer coefficientssoilThermalTransfer
        soilTransfer transfer;
        transfer.d.resize(depth);
        transfer.k.resize(depth);
        
        // loop through each depth
        for (int d=0; d<depth; ++d) {
            
            transfer.d[d] = (b[d]*K_nsat[d]*std::abs(psi_nsat[d])/soil_q[d])
                            *std::pow(soil_q[d]/porosity[d],(b[d]+3.));
            transfer.k[d] = K_nsat[d]*std::pow(soil_q[d]/porosity[d],(2.*b[d]+3.));
        }
                
        return transfer;
    }
    
    // set soil type properties at each depth
    // soil type from USDA 11-category + peat
    // 01 = sand
    // 02 = loamy sand
    // 03 = sandy loam
    // 04 = silty loam
    // 05 = loam
    // 06 = sandy clay loam
    // 07 = silty clay loam
    // 08 = clay loam
    // 09 = sandy clay
    // 10 = silty clay
    // 11 = clay
    // 12 = peat
    soilProperties soilTypeProperties(const std::vector<int>& soil_type, const int depth, const int src) {
        
        
        std::vector<double>b_list;   // exponent (unitless)
        std::vector<double>psi_list; // saturation moisture potential (m)
        std::vector<double>por_list; // porosity (volume/volume)
        std::vector<double>res_list; // residual moisture (volume/volume)
        std::vector<double>K_list;   // hydraulic conductivity (m/s)
        std::vector<double>Ci_list;  // volumetric heat capacity (J/m^3/K)
        
        // struct to hold soil properties
        soilProperties properties;
        properties.b.resize(depth);
        properties.psi_sat.resize(depth);
        properties.porosity.resize(depth);
        properties.K_sat.resize(depth);
        properties.Ci.resize(depth);
        
        switch (src) {
            
            // Clapp and Hornberger (1974)
            case 1:
                b_list   = { 4.05, 4.38,  4.90,  5.30,  5.39,  7.12,
                             7.75, 8.52, 10.40, 10.40, 11.40,  7.75};
                
                psi_list = { -0.121, -0.090, -0.218, -0.786, -0.478, -0.299,
                             -0.356, -0.630, -0.153, -0.490, -0.405, -0.356};
                
                por_list = { 0.395, 0.410, 0.435, 0.485, 0.451, 0.420,
                             0.477, 0.476, 0.426, 0.492, 0.482, 0.863};
                
                res_list = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
                
                K_list   = { 1.76e-4, 1.56e-4, 3.41e-5, 7.20e-6, 7.00e-6, 6.30e-6,
                             1.70e-6, 2.50e-6, 2.20e-6, 1.00e-6, 1.30e-6, 8.00e-6};
                
                Ci_list  = { 1.47e6, 1.41e6, 1.34e6, 1.27e6, 1.21e6, 1.18e6,
                             1.32e6, 1.23e6, 1.18e6, 1.15e6, 1.09e6, 0.84e6};
            
            // Cosby et al. (1984)
            case 2:
                b_list   = { 2.79, 4.26,  4.74,  5.33,  5.25, 6.77,
                             8.72, 8.17, 10.73, 10.39, 10.55, 7.75};
                
                psi_list = {-0.023, -0.018, -0.032, -0.066, -0.047, -0.031,
                            -0.060, -0.041, -0.027, -0.045, -0.053, -0.356};
                
                por_list = {0.339, 0.421, 0.434, 0.476, 0.439, 0.404,
                            0.464, 0.465, 0.406, 0.468, 0.468, 0.863};
                
                res_list = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
                
                K_list   = {1.60e-5, 9.52e-6, 6.19e-6, 4.73e-6, 5.12e-6, 5.78e-6,
                            4.11e-6, 4.45e-6, 7.12e-5, 3.43e-6, 2.99e-6, 8.00e-6};
                
                Ci_list  = {1.47e6, 1.41e6, 1.34e6, 1.27e6, 1.21e6, 1.18e6,
                            1.32e6, 1.23e6, 1.18e6, 1.15e6, 1.09e6, 0.84e6};
            
            // Rawls and Brakensiek (1982) [uses geomtric means for psi and b
            case 3:
                b_list   = { 1.69, 2.11, 3.11, 4.55, 4.74, 4.00,
                             5.15, 6.62, 5.95, 7.87, 7.63, 7.75};
                
                psi_list = { -0.073, -0.087, -0.147, -0.112, -0.208, -0.281,
                             -0.259, -0.326, -0.292, -0.342, -0.373, -0.356};
                
                por_list = { 0.437, 0.437, 0.453, 0.463, 0.501, 0.398,
                             0.464, 0.471, 0.430, 0.479, 0.475, 0.863};
                
                res_list = { 0.020, 0.035, 0.041, 0.027, 0.015, 0.068,
                             0.075, 0.040, 0.109, 0.056, 0.090, 0.1763};
                
                K_list   = { 5.83E-5, 1.70e-5, 7.19e-6, 1.89e-6, 3.67e-6, 1.19e-6,
                             6.39e-7, 4.17e-7, 3.33e-7, 2.50e-7, 1.67e-7, 8.00e-6};
                
                Ci_list  = { 1.47e6, 1.41e6, 1.34e6, 1.27e6, 1.21e6, 1.18e6,
                             1.32e6, 1.23e6, 1.18e6, 1.15e6, 1.09e6, 0.84e6};
        }
         
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
