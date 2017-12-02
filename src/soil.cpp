//
//  soil.hpp
//  
//  This namespace handles functions related to soil 
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include "soil.hpp"
#include "constants.hpp"
#include <cmath> 
#include <tuple>

namespace soil {
    
    // compute surface mixing ratio
    double surfaceMixingRatio(const double psi_nsat, const double porosity, 
                              const double b, const double sfc_T, 
                              const double sfc_q, const double atm_p) {
        
        double psi_n = psi_nsat*std::pow((porosity/sfc_q),b);
        double h = std::exp(Constants::grav*psi_n/(Constants::Rv*sfc_T));
        double e = 6.1078*std::exp(17.269*(sfc_T-273.15)/(sfc_T-35.86));
        double satHum = 0.622*(e/(atm_p-0.378*e));
        double specHum = h*satHum;
        
        return specHum/(1-specHum);
    }
    
    // compute average soil conductivity/diffusivity
    double soilThermalTransfer(const double *psi_nsat, const double *porosity, 
                               const double *soil_q, const double *b, 
                               const double *Ci, const int depth, const int flag) {
        
        double transfer = 0;
        double psi_n, pf;
        double temp;
        
        // loop through each depth
        for (int d=0; d<depth; ++d) {
            psi_n = 100.*psi_nsat[d]*std::pow((porosity[d]/soil_q[d]),b[d]);
            pf = std::log10(std::abs(psi_n));
            if (pf <= 5.1) {
                transfer += 418.46*std::exp(-(pf+2.7));
            } else {
                transfer += 0.172;
            }
                        
            // convert to diffusivity if flag==1
            if (flag==1) {
                double heatCap = (1-porosity[d])*Ci[d] + soil_q[d]*Constants::Ci_wat;
                transfer = transfer / heatCap;
            }
        }
        
        return transfer/depth;
    }
    
    // compute average soil moisture transfer
    std::tuple<double, double> soilMoistureTransfer(const double *psi_nsat, const double *K_nsat, 
                                                    const double *porosity, const double *soil_q, 
                                                    const double *b, const int depth) {
        
        double transfer_h=0, transfer_d=0;
        
        // loop through each depth
        for (int d=0; d<depth; ++d) {
            transfer_d += -(b[d]*K_nsat[d]*psi_nsat[d]/soil_q[d])*std::pow(soil_q[d]/porosity[d],b[d]+3.);
            transfer_h += K_nsat[d]*std::pow(soil_q[d]/porosity[d],2*b[d]+3.);
        }
        
        return std::make_tuple(transfer_d/depth, transfer_h/depth); 
    }
    
};