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

/**
 * This is the main UtahLSM class.
 */

class UtahLSM {
    
    private:
        
        // input variables
        bool first;
        double dt, z_o, z_t, z_m, z_s, atm_p, atm_ws, atm_T, atm_q;
        int nsoilz;
        std::vector<double> &soil_z;
        std::vector<int> &soil_type;
        std::vector<double> &soil_T;
        std::vector<double> &soil_T_last;
        std::vector<double> &soil_q;
        std::vector<double> &soil_q_last;
        int julian_day;
        double utc, latitude, longitude, albedo, emissivity, R_net;
        int comp_rad;
        double &zeta_m, &zeta_s,  &zeta_o,  &zeta_t;
        double &ustar,  &flux_wT, &flux_wq, &flux_gr;
        
        // local variables
        double surf_T_last;
        double surf_q_last;

        std::vector<double> b;
        std::vector<double> psi_sat;
        std::vector<double> porosity;
		std::vector<double> residual;
		std::vector<double> K_sat;
		std::vector<double> Ci;
        
        // functions
        void setSoilProperties();
        void computeFluxes(double,double);
        void computeRadiation();
        void solveSEB();
        void solveSMB();
        void solveMoisture();
        void solveDiffusion(int);
        double computeSEB(double);
        double computeDSEB(double);
        double computeSMB(double);
        double computeDSMB(double);
        
    public :
        UtahLSM(bool, double, double, double, double, double,
                double, double, double, double,
                int, std::vector<double>&, std::vector<int>&, 
                std::vector<double>&, std::vector<double>&,
                std::vector<double>&, std::vector<double>&, 
                int, double, double, double,
                double, double, double,int,
                double&,double&,double&,double&,
                double&,double&,double&,double&);                
};

#endif
