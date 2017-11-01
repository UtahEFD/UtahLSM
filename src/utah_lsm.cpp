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
    
    // set constants
    vonk = 0.4;               // von Karman constant
    grav = 9.81;              // gravity
    pi   = 3.14159265358979;  // pi
    sb   = 5.6697e-8;         // Stefan-Boltzmann constant
    sc   = 1.127;             // solar constant [K-m/s]
    rhoA = 1.204;             // density of air [kg/m^3]
    rhoW = 1000.0;            // density of water [kg/m^3]
    Rv   = 461.4;             // gas constant for water vapor [J/kg-K]
    Lv   = 2.45e6;            // latent heat of vaporization [J/kg]
}

// Initialize UtahLSM elements using input data
void UtahLSM :: init() {
    
    std::cout<<"Initializing UtahLSM";
    std::cout<<"##############################################################"<<std::endl;
}