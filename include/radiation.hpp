//
//  radiation.hpp
//  
//  This class computes surface radiation
//
//  Created by Jeremy Gibbs on 11/6/17.
//

#ifndef RADIATION_HPP
#define RADIATION_HPP

/**
 * This class computes the surface radiation.
 */

class Radiation {
    
    public:
        
        Radiation(const double,const double,const double,const double);
        
        double computeNet(const double,const double,const double);
        
    private:
        
        // site properties
        double latitude, longitude, albedo, emissivity;
        
        // functions for each component
        double longwaveIn();
        double longwaveOut(const double,const double);
        double shortwaveIn(const double,const double,const double,const double);
        double shortwaveOut(const double,const double);
    
};

#endif
