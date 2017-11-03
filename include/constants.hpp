//
//  constants.hpp
//  
//  This class stores frequently used constants
//
//  Created by Jeremy Gibbs on 11/3/17.
//

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace Constants {
    const double vonk = 0.4;               // von Karman constant
    const double grav = 9.81;              // gravity
    const double pi   = 3.14159265358979;  // pi
    const double sb   = 5.6697e-8;         // Stefan-Boltzmann constant
    const double sc   = 1.127;             // solar constant [K-m/s]
    const double rhoA = 1.204;             // density of air [kg/m^3]
    const double rhoW = 1000.0;            // density of water [kg/m^3]
    const double Rv   = 461.4;             // gas constant for water vapor [J/kg-K]
    const double Lv   = 2.45e6;            // latent heat of vaporization [J/kg]
};

#endif