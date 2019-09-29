/*
 * UtahLSM
 * 
 * Copyright (c) 2019 Jeremy A. Gibbs
 * Copyright (c) 2019 Pete Willemsen
 * Copyright (c) 2019 Rob Stoll
 * Copyright (c) 2019 Eric Pardyjak
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
                     const double albedo,const double emissivity,
                     const int nx,const int ny) : 
                     latitude(latitude), longitude(longitude), 
                     albedo(albedo), emissivity(emissivity),
                     nx(nx), ny(ny) {}

// Computes the surface net radiation
std::vector<double> Radiation::computeNet(const double julian_day, const double time_utc, 
                                          const std::vector<double> sfc_T) {
    
    // Compute incoming shortwave radiation
    std::vector<double> shortwave_in = shortwaveIn(julian_day, time_utc, latitude, longitude);
    
    // Compute outgoing shortwave radiation
    std::vector<double> shortwave_out = shortwaveOut(albedo, shortwave_in);
    
    // Compute outgoing longwave radiation
    std::vector<double> longwave_out = longwaveOut(emissivity, sfc_T);

    // Compute net radiation
    std::vector<double> net_r(ny*nx,0.0);
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            int id = i + j*nx;
            
            // Compute incoming longwave radiation (current hack is net longwave of -50)
            double longwave_in = longwave_out[id] - 50.0;

            net_r[id] = shortwave_in[id] - shortwave_out[id] + longwave_in - longwave_out[id];
        }
    }
    return net_r;
}

// Compute outgoing longwave radiation
std::vector<double> Radiation::longwaveOut(const double emissivity,const std::vector<double> sfc_T) {
    
    std::vector<double> longwave_out(ny*nx,0.0);
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            int id = i + j*nx;
            longwave_out[id] = emissivity * c::sb * std::pow(sfc_T[id],4.);
        }
    }
    return longwave_out;
}

// Compute incoming shortwave radiation
std::vector<double> Radiation::shortwaveIn(const double julian_day,const double time_utc,
                                           const double latitude,const double longitude) {
    
    // local variables
    std::vector<double> shortwave_in(ny*nx,0.0);
    double declination   = 23.45*(c::pi/180.0)*std::cos(2.0*c::pi*(julian_day-173)/365.25);
    double sin_elevation = std::sin(latitude)*std::sin(declination) - 
                           std::cos(latitude)*std::cos(declination) * 
                           std::cos( (2*c::pi*time_utc/(24.0*3600.0)) - longitude );
                           
    if (sin_elevation > 0) {
        double transmissivity = (0.6 + 0.2*sin_elevation);
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                int id = i + j*nx;
                shortwave_in[id] = c::sc * transmissivity * sin_elevation;
            }
        }
    }
    
    return shortwave_in;
}

// Compute outgoing shortwave radiation
std::vector<double> Radiation::shortwaveOut(const double albedo, const std::vector<double> shortwave_in) {
    
    // make a copy of shortwave in
    std::vector<double> shortwave_out = shortwave_in;
    
    // multiply shortwave in by albedo
    std::transform(shortwave_out.begin(), shortwave_out.end(), shortwave_out.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, albedo));
    
    return shortwave_out;
}