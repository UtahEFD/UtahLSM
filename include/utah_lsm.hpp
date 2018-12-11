//
//  utah_lsm.hpp
//  
//  This class handles the Utah LSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#ifndef UTAHLSM_HPP
#define UTAHLSM_HPP

#include "input.hpp"
#include <string>
#include <vector>

/**
 * This is the main UtahLSM class.
 */

class Input;

class UtahLSM {
    
    private:
        
        // flux pointers
        double &ustar, &flux_wT, &flux_wq;

        // namelist length section
        double z_o, z_t, z_m, z_s;
        
        // namelist soil section
        int soil_param, soil_model, nsoilz;
        std::vector<int> soil_type;
        std::vector<double> soil_z;
        std::vector<double> soil_T;
        std::vector<double> soil_q;
    
        // local variables
        double flux_gr;
        double zeta_m, zeta_s, zeta_o, zeta_t;
        double surf_T_last, surf_q_last; 
        std::vector<double> soil_T_last;
        std::vector<double> soil_q_last;
    
        // namelist radiation section
        int utc_start, julian_day, comp_rad;
        double albedo, emissivity, latitude, longitude;

        // soil properties
        std::vector<double> b;
        std::vector<double> psi_sat;
        std::vector<double> porosity;
        std::vector<double> residual;
        std::vector<double> K_sat;
        std::vector<double> Ci;

        // atmospheric data
        double atm_U, atm_T, atm_q, atm_p, R_net;
    
        // time data
        bool first=true;
        double tstep, runtime, utc=0;
    
        // internal functions
        void setSoilProperties();
        void computeRadiation();
        void computeFluxes(double,double);
        void solveSEB();
        void solveSMB();
        void solveDiffusion(int);
        double computeSEB(double);
        double computeDSEB(double);
    
    public :
        
        UtahLSM(Input*,double&,double&,double&);
    
        // external functions
        void updateFields(double,double,double,double,double,double);
        void run();
};

#endif
