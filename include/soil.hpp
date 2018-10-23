//
//  soil.hpp
//  
//  This namespace handles functions related to soil
//
//  Created by Jeremy Gibbs on 11/6/17.
//

#ifndef SOIL_HPP
#define SOIL_HPP

#include <tuple>
#include <vector>

/**
 * This class handles tasks associated with soil.
 */

namespace soil {
	
	struct moistureTransfer {
        std::vector<double> d; 
        std::vector<double> k;
    };
    
    struct thermalTransfer {
        std::vector<double> d;
        std::vector<double> k;
    };
    
    struct properties {
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
    
    double surfaceWaterContent(const double,const double,const double,
                               const double,const double,const double,
                               const double,const int);
    
    double heatCapacity(const double, const double, const double);
    
    double waterPotential(const double, const double, const double,
                          const double, const double, const int);
    
    std::vector<double> waterPotential(const std::vector<double>&, const std::vector<double>&,
                                       const std::vector<double>&, const std::vector<double>&,
                                       const std::vector<double>&, const int, const int);
    
    thermalTransfer thermalTransfer(const std::vector<double>&, const std::vector<double>&,
                                        const std::vector<double>&, const std::vector<double>&,
                                        const std::vector<double>&, const int);
    
    moistureTransfer moistureTransfer(const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const std::vector<double>&,
                                      const int, const int);
    
    properties properties(const std::vector<int>&, const int, const int);
};

#endif
