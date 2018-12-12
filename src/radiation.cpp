//
//  radiation.cpp
//  
//  This namespace handles functions related to radiation 
//
//  Created by Jeremy Gibbs on 10/31/17.
//
#include <cmath> 
#include "radiation.hpp"
#include "constants.hpp"

namespace {
    namespace c = Constants;
}

namespace radiation {
    
    // compute incoming longwave radiation
    double longwaveIn() {
        
        return 0;
    }
    
    // compute outgoing longwave radiation
    double longwaveOut(const double emissivity, const double sfc_T) {
        
        return emissivity * c::sb * std::pow(sfc_T,4.);
    }
    
    // compute incoming shortwave radiation
    double shortwaveIn(const double julian_day,const double time_utc,
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
    
    // compute outgoing shortwave radiation
    double shortwaveOut(const double albedo, const double shortwave_in) {
        
        return albedo * shortwave_in;
    }
    
    // compute net radiation
    double net(const double sfc_T,const double emissivity, 
               const double julian_day,const double time_utc,
               const double latitude,const double longitude,
               const double albedo) {
        
        double shortwave_in  = shortwaveIn(julian_day, time_utc,latitude, longitude);
        double shortwave_out = shortwaveOut(albedo,shortwave_in);
        double longwave_out  = longwaveOut(emissivity,sfc_T);
        
        // current hack to estimate net longwave of -50.
        double longwave_in   = longwave_out - 50.0;
        
        return shortwave_in - shortwave_out + longwave_in - longwave_out;
    }
    
};
