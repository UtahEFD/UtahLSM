//
//  utah_lsm.hpp
//  
//  This class handles the Utah LSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#ifndef UTAHLSM_HPP
#define UTAHLSM_HPP

#include <string>
//#include <netcdf>
//#include "input.hpp"
//using namespace netCDF;
//using namespace netCDF::exceptions;

class UtahLSM {
    
    private:
        float vonk;
        float grav;
        float pi;       
        
    public :
        UtahLSM();
                
        // Functions
        void init();
};

#endif