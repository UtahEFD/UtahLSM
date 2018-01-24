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

namespace soil {
    
    double surfaceMixingRatio(const double,const double, const double,
                              const double,const double, const double);
                              
    std::vector<double> soilThermalTransfer(const std::vector<double>&, const std::vector<double>&, 
                                            const std::vector<double>&, const std::vector<double>&, 
                                            const std::vector<double>&, const int, const int);
    
    struct soilTransfer {
        std::vector<double> transfer_d;
        std::vector<double> transfer_h;
    };
    
    soilTransfer soilMoistureTransfer(const std::vector<double>&, const std::vector<double>&, 
                                      const std::vector<double>&, const std::vector<double>&,
                                      const std::vector<double>&, const int);
};

#endif