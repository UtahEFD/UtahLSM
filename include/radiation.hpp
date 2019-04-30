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
