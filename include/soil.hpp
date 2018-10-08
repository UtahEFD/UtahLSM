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
	
	struct soilTransfer {
        std::vector<double> d; 
        std::vector<double> k;
    };
    
    struct soilProperties {
	    
	    std::vector<double> b;
	    std::vector<double> psi_sat;
	    std::vector<double> porosity;
	    std::vector<double> K_sat;
	    std::vector<double> Ci;
    };
    
    double surfaceMixingRatio(const double,const double, const double,
                              const double,const double, const double);
    
    double surfaceSoilMoisture(const double,const double, const double,
                               const double,const double, const double);
                              
    std::vector<double> soilThermalTransfer(const std::vector<double>&, const std::vector<double>&, 
                                            const std::vector<double>&, const std::vector<double>&, 
                                            const std::vector<double>&, const int, const int);
    
    soilTransfer soilMoistureTransfer(const std::vector<double>&, const std::vector<double>&, 
                                      const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const int);
    
    soilProperties soilTypeProperties(const std::vector<int>&, const int);
};

#endif
