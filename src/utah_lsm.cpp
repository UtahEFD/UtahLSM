//
//  utah_lsm.cpp
//  
//  This class handles the Utah LSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include "utah_lsm.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <mpi.h>
#include <fstream>
#include <stdexcept>
//#include <netcdf>
//using namespace netCDF;
//using namespace netCDF::exceptions;

UtahLSM::UtahLSM() {
}

// Initialize UtahLSM elements using input data
void UtahLSM :: init() {
    
    std::cout<<"Initializing UtahLSM";
    std::cout<<"##############################################################"<<std::endl;
}