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

class UtahLSM {
    
    private:
    
        // local variables
        double lastFlux, lastTemp;
        double tempConverged, fluxConverged;
        
        // input variables
        int nsoilz, julian_day;
        double z_o, z_t, z_m, z_s, atm_ws, atm_T, atm_q, atm_p;
        double utc, latitude, longitude, albedo, emissivity, R_net;
        double *soil_z, *soil_T, *soil_q;
        double *porosity, *psi_nsat, *K_nsat, *b, *Ci;
        double &zeta_m,&zeta_s,&zeta_o,&zeta_t;
        double &ustar, &flux_wT, &flux_wq;
    
        // functions
        void computeFluxes();
        void solveSEB();
        bool convergeTemp();
        bool convergeFlux();   
        
    public :
        UtahLSM(double, double, double, double,
                double, double, double, double,
                int, double*, double*, double*,
                double*, double*, double*, double*, double*,
                int, double, double, double,
                double, double, double,
                double&,double&,double&,double&,
                double&,double&,double&);                
};

#endif