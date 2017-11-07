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
        // variables
        int nsoilz, julian_day;
        double z_o, z_t, z_m, z_s, atm_ws, atm_T, atm_q, atm_p, L;
        double utc, latitude, longitude, albedo, emissivity, R_net;
        double *soil_z, *soil_T, *soil_q;
        double *porosity, *psi_nsat, *K_nsat, *b, *Ci;
        double &phiM, &psiM, &psiM0, &phiH, &psiH, &psiH0;
        double &ustar, &flux_wT, &flux_wq;
    
        // functions
        void computeFluxes();       
        
    public :
        UtahLSM(double, double, double, double,
                double, double, double, double,
                int, double*, double*, double*,
                double*, double*, double*, double*, double*,
                int, double, double, double,
                double, double, double,
                double&,double&,double&,
                double&,double&,double&,
                double&,double&,double&);                
};

#endif