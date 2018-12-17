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
#include "output.hpp"
#include <string>
#include <vector>

/**
 * This is the main UtahLSM class.
 */

class Input;
class Output;
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
        double flux_gr=0, L=0;
        double zeta_m=0, zeta_s=0, zeta_o=0, zeta_t=0;
        double surf_T_last=0, surf_q_last=0;
        std::vector<double> soil_T_last;
        std::vector<double> soil_q_last;
    
        // namelist radiation section
        int utc_start, julian_day, comp_rad;
        double albedo, emissivity, latitude, longitude;
    
        // namelist output section
        int saveOutput = false;
        std::vector<std::string> outputFields;

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
        double tstep=0, runtime=0, utc=0;
    
        // local output information
        Output* output;
    
        int output_counter=0;
    
        struct attributes1D {
            double data;
            std::string name;
            std::string long_name;
            std::string units;
        };
    
        struct attributes2D {
            std::vector<double> data;
            std::string name;
            std::string long_name;
            std::string units;
        };
    
        std::map<std::string,attributes1D> map1D;
        std::map<std::string,attributes2D> map2D;
    
        std::vector<attributes1D> fieldsToSave1D;
        std::vector<attributes2D> fieldsToSave2D;
    
        std::vector<attributes1D> attList1D;
        std::vector<attributes2D> attList2D;
    
        std::vector<NcDim> dim_1D;
        std::vector<NcDim> dim_2D;
    
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
        void save();
};

#endif
