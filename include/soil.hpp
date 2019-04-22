//
//  soil.hpp
//  
//  This namespace handles functions related to soil
//
//  Created by Jeremy Gibbs on 11/6/17.
//

#ifndef SOIL_HPP
#define SOIL_HPP

#include <vector>
include <stdexcept>

/**
 * This class handles tasks associated with soil.
 */

namespace soil {
	
	struct MoistureTransfer {
        std::vector<double> d; 
        std::vector<double> k;
    };
    
    struct ThermalTransfer {
        std::vector<double> d;
        std::vector<double> k;
    };
    
    struct Properties {
	    std::vector<double> b;
	    std::vector<double> psi_sat;
	    std::vector<double> porosity;
        std::vector<double> residual;
	    std::vector<double> K_sat;
	    std::vector<double> Ci;
    };
    
    double heatCapacity(const double, const double, const double);
    
    double surfaceMixingRatio(const double,const double, const double,
                              const double,const double, const double,
                              const double, const int);
    
    double surfaceWaterContent(const double, const double, const double,
                               const double, const double, const int);
    
    double surfaceWaterContentEstimate(const double,const double,const double,
                                       const double,const double,const double,
                                       const double,const int);
    
    double heatCapacity(const double, const double, const double);
    
    double waterPotential(const double, const double, const double,
                          const double, const double, const int);
    
    std::vector<double> waterPotential(const std::vector<double>&, const std::vector<double>&,
                                       const std::vector<double>&, const std::vector<double>&,
                                       const std::vector<double>&, const int, const int);
    
    ThermalTransfer thermalTransfer(const std::vector<double>&, const std::vector<double>&,
                                    const std::vector<double>&, const std::vector<double>&,
                                    const std::vector<double>&, const std::vector<double>&,
                                    const int, const int);
    
    MoistureTransfer moistureTransfer(const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const std::vector<double>&,
                                      const int, const int);
    
    Properties properties(const std::vector<int>&, const int, const int);
};

#endif
