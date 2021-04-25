/*
 * UtahLSM
 * 
 * Copyright (c) 2021 Jeremy A. Gibbs
 * Copyright (c) 2021 Rob Stoll
 * Copyright (c) 2021 Eric Pardyjak
 * Copyright (c) 2021 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "radiation.hpp"

#include <cmath> 

#include "constants.hpp"

namespace {
    namespace c = constants;
}

// Constructor for Radiation class
Radiation::Radiation(const double latitude,const double longitude,
                     const double albedo,const double emissivity) : 
                     latitude(latitude), longitude(longitude), 
                     albedo(albedo), emissivity(emissivity) {}

// Computes the surface net radiation
double Radiation::computeNet(const double julian_day, const double time_utc, const double sfc_T) {
    
    // Compute incoming shortwave radiation
    double shortwave_in = shortwaveIn(julian_day, time_utc, latitude, longitude);
    
    // Compute outgoing shortwave radiation
    double shortwave_out = shortwaveOut(albedo, shortwave_in);
    
    // Compute outgoing longwave radiation
    double longwave_out = longwaveOut(emissivity, sfc_T);
    
    // Compute incoming longwave radiation (current hack is net longwave of -50)
    double longwave_in = longwave_out - 50.0;
    
    // Compute net radiation
    return shortwave_in - shortwave_out + longwave_in - longwave_out;
}

// Compute incoming longwave radiation
double Radiation::longwaveIn() {
    
    return 0;
}

// Compute outgoing longwave radiation
double Radiation::longwaveOut(const double emissivity, const double sfc_T) {
    
    return emissivity * c::sb * std::pow(sfc_T,4.);
}

// Compute incoming shortwave radiation
double Radiation::shortwaveIn(const double julian_day,const double time_utc,
                              const double latitude,const double longitude) {
    
    // local variables
    double shortwave_in = 0;
    double declination = 23.45*(c::pi/180.0)*std::cos(2.0*c::pi*(julian_day-173)/365.25);
    double sin_elevation = std::sin(latitude)*std::sin(declination) - 
                           std::cos(latitude)*std::cos(declination) * 
                           std::cos( (2*c::pi*time_utc/(24.0*3600.0)) - longitude );
                           
    if (sin_elevation > 0) {
        double transmissivity = (0.6 + 0.2*sin_elevation);
        shortwave_in = c::sc * transmissivity * sin_elevation;
    }
    
    return shortwave_in;
}

// Compute outgoing shortwave radiation
double Radiation::shortwaveOut(const double albedo, const double shortwave_in) {
    
    return albedo * shortwave_in;
}