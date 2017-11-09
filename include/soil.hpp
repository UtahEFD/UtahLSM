//
//  soil.hpp
//  
//  This namespace handles functions related to soil
//
//  Created by Jeremy Gibbs on 11/6/17.
//

#ifndef SOIL_HPP
#define SOIL_HPP

namespace soil {
    
    double surfaceMixingRatio(const double,const double, const double,
                              const double,const double, const double);
                              
    double soilThermalTransfer(const double*, const double*, const double*,
                               const double*, const double*, const int, const int);
};

#endif