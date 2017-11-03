//
//  utah_lsm.cpp
//  
//  This class handles the Utah LSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include "utah_lsm.hpp"
#include "constants.hpp"
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

UtahLSM::UtahLSM(double z_o,double z_t,double z_m,double z_s, 
                 double air_u,double air_T,double air_q, 
                 int nsoilz,double* soil_z,double* soil_T,double* soil_q,
                 double* porosity,double* psi_nsat,double* K_nsat,double* b,double* Ci,
                 int julian_day, double utc, double latitude, double longitude,
                 double albedo, double emissivity, double R_net,
                 double* phi, double* psi, double* psi0,
                 double* phiH, double* psiH, double* psiH0,
                 double* ustar, double* flux_wT, double* flux_wq) {
    
        
}