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
    
    // compute heat capacity
    double heatCapacity(const double porosity, const double Ci, const double soil_q) {
        return (1-porosity)*Ci + soil_q*c::Ci_wat + (porosity-soil_q)*c::Cp_air;
    }
    
    // compute surface mixing ratio
    double surfaceMixingRatio(const double psi_sat, const double porosity,
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
    
    // compute soil surface moisture
    double surfaceWaterContent(const double psi, const double psi_sat,
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
    
    // estimate soil surface moisture from surface mixing ratio
    double surfaceWaterContentEstimate(const double psi_sat, const double porosity,
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
    
    // compute soil thermal conductivity/diffusivity
    struct ThermalTransfer thermalTransfer(const std::vector<double> &psi_sat,
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
    
    // compute average soil moisture transfer
    struct MoistureTransfer moistureTransfer(const std::vector<double>& psi_sat,
                                             const std::vector<double>& K_sat,
                                             const std::vector<double>& porosity,
                                             const std::vector<double>& residual,
                                             const std::vector<double>& soil_q,
                                             const std::vector<double>& b,
                                             const int depth, const int model) {
        
        // local variables
        double Se;
        struct MoistureTransfer transfer;
        transfer.k.resize(depth);
        transfer.d.resize(depth);

        // loop through each depth
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
    
    // compute soil water potential (single level)
    double waterPotential(const double psi_sat, const double porosity,
                          const double residual, const double soil_q,
                          const double b, const int model) {
        
        // local variables
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
    
    // compute soil water potential (full column)
    std::vector<double> waterPotential(const std::vector<double>& psi_sat,
                                       const std::vector<double>& porosity,
                                       const std::vector<double>& residual,
                                       const std::vector<double>& soil_q,
                                       const std::vector<double>& b,
                                       const int depth, const int model) {
        
        // local variables
        double Se;
        std::vector<double> psi(depth);
        
        // loop through each depth
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
    struct Properties properties(const std::vector<int>& soil_type, const int depth, const int src) {
        
        std::vector<double>b_list;   // exponent (unitless)
        std::vector<double>psi_list; // saturation moisture potential (m)
        std::vector<double>por_list; // porosity (volume/volume)
        std::vector<double>res_list; // residual moisture (volume/volume)
        std::vector<double>K_list;   // hydraulic conductivity (m/s)
        std::vector<double>Ci_list;  // volumetric heat capacity (J/m^3/K)
        
        // struct to hold soil properties
        struct Properties properties;
        properties.b.resize(depth);
        properties.psi_sat.resize(depth);
        properties.porosity.resize(depth);
        properties.residual.resize(depth);
        properties.K_sat.resize(depth);
        properties.Ci.resize(depth);
        
        switch (src) {
            
            // Clapp and Hornberger (1974)
            case 1:
            {
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
                break;
            }
            // Cosby et al. (1984)
            case 2:
            {
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
                break;
            }
            // Rawls and Brakensiek (1982) [uses arithmetic means for psi and b]
            case 3:
            {
                b_list   = { 1.44, 1,81, 2.65, 3.97, 4.27, 3.13,
                             4.13, 5.65, 4.48, 6.67, 6.06, 7.75};
                
                psi_list = { -0.160, -0.206, -0.302, -0.401, -0.509, -0.594,
                             -0.564, -0.703, -0.795, -0.765, -0.856, -0.356};
                
                por_list = { 0.437, 0.437, 0.453, 0.463, 0.501, 0.398,
                             0.464, 0.471, 0.430, 0.479, 0.475, 0.863};
                
                res_list = { 0.020, 0.035, 0.041, 0.027, 0.015, 0.068,
                             0.075, 0.040, 0.109, 0.056, 0.090, 0.1763};
                
                K_list   = { 5.83E-5, 1.70e-5, 7.19e-6, 1.89e-6, 3.67e-6, 1.19e-6,
                             6.39e-7, 4.17e-7, 3.33e-7, 2.50e-7, 1.67e-7, 8.00e-6};
                
                Ci_list  = { 1.47e6, 1.41e6, 1.34e6, 1.27e6, 1.21e6, 1.18e6,
                             1.32e6, 1.23e6, 1.18e6, 1.15e6, 1.09e6, 0.84e6};
                break;
            }
            // Wosten (Cabauw specific)
            case 4:
            {
                b_list   = { 9.09, 6.25, 2.63};
                
                psi_list = { -0.513, -1.053, -0.971};
                
                por_list = { 0.590, 0.560, 0.890};
                
                res_list = { 0.010, 0.010, 0.0};
                
                K_list   = { 5.2e-7, 1.2e-7, 1.2e-7};
                
                Ci_list  = { 1.2e6, 1.2e6, 1.2e6};
                
                break;
            }
            default:
            {
                throw std::invalid_argument("soil_param must be set to 1, 2, or 3");
                break;
            }
        }
         
        // loop through each depth to assign soil type properties
        int soil;
        for (int d=0; d<depth; ++d) {
	        
	        soil = soil_type[d] - 1;
	        
            properties.b[d]        = b_list[soil];
            properties.psi_sat[d]  = psi_list[soil];
            properties.porosity[d] = por_list[soil];
            properties.residual[d] = res_list[soil];
            properties.K_sat[d]    = K_list[soil];
            properties.Ci[d]       = Ci_list[soil];
        }
        return properties;
    }  
};
