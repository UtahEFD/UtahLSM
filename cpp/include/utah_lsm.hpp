/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2023 Jeremy A. Gibbs
 * Copyright (c) 2017–2023 Rob Stoll
 * Copyright (c) 2017–2023 Eric Pardyjak
 * Copyright (c) 2017–2023 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#ifndef UTAHLSM_HPP
#define UTAHLSM_HPP

#include <netcdf>
#include <string>
#include <vector>

using namespace netCDF;
using namespace netCDF::exceptions;

class Input;
class Logger;
class Output;
class Settings;
class Soil;
class Surface;
class Radiation;

/**
 * This is the main UtahLSM class.
 * 
 * This class is responsible for solving the surface energy 
 * and moisture budgets to compute the associated fluxes 
 * of momentum, heat, and moisture.
 */
class UtahLSM {
    
    public:
        
        /**
         * Constructs a UtahLSM object.
         *
         * @param[in]     settings Settings object
         * @param[in]     input Input object
         * @param[in]     output Output object
         * @param[in,out] ustar friction velocity
         * @param[in,out] flux_wT kinematic heat flux
         * @param[in,out] flux_wq kinematic moisture flux
         * @param[in]     j y-index of UtahLSM column if in a grid
         * @param[in]     i x-index of UtahLSM column if in a grid
         */
        UtahLSM(Settings* settings, Input* input, Output* output, double& ustar, 
                double& flux_wT,double& flux_wq, int j=0, int i=0);
    
        /**
         * Update atmospheric quantities prior to solving.
         *
         * @param[in] dt current time-step
         * @param[in] u wind speed
         * @param[in] T temperature
         * @param[in] q mixing ratio
         * @param[in] p pressure
         * @param[in] rad optional net radiation
         */
        void updateFields(double dt,double u,double T,double q,double p,double rad);
        
        /**
         * Run UtahLSM.
         */
        void run();

        /**
         * Save current fields to the output file.
         *
         * @param[in] output Output object
         */
        void save(Output* output);
    
    private:

        // Classes 
        Soil* soil;
        Surface* sfc;
        Radiation* radiation;
        Logger* logger;

        // Diffusion time step restrictions
        double dt_T; ///< diffusion time step for temperature
        double dt_q; ///< diffusion time step for moisture

        // Fields to compute
        double &ustar;   ///< reference to friction velocity
        double &flux_wT; ///< reference to kinematic heat flux
        double &flux_wq; ///< reference to kinematic moisture flux
        double flux_sh;  ///< sensible heat flux
        double flux_lh;  ///< latent heat flux

        // Input time section
        int step_seb;   ///< steps between calls to energy balance
        int step_dif;   ///< steps between calls to diffusion
        int utc;        ///< time in UTC at simulation start
        int julian_day; ///< Julian day at simulation start

        // Input grid section
        int nx; ///< number of points in x-direction
        int ny; ///< number of points in y-direction
        
        // Input surface section
        double z_o;        ///< roughness height for momentum
        double z_t;        ///< roughness height for scalars
        double z_m;        ///< measurement height for momentum
        double z_s;        ///< measurement height for scalars
        double albedo;     ///< surface albedo at site
        double emissivity; ///< surface emissivity at site
        
        // Input soil section
        int soil_param;             ///< soil parameter set
        int soil_model;             ///< soil model
        int nsoilz;                 ///< number of soil levels
        std::vector<int> soil_type; ///< soil type at each level
        std::vector<double> soil_z; ///< soil depth at each level
        std::vector<double> soil_T; ///< soil temperature at each level
        std::vector<double> soil_q; ///< soil moisture at each level

        // Input radiation section
        int comp_rad;      ///< flag whether to run radiation model or not
        double latitude;   ///< latitude at site
        double longitude;  ///< longitude at site
    
        // Input output section
        int save_output = false;                ///< flag whether to save output
        std::vector<std::string> output_fields; ///< list of fields to save

        // Local surface variables
        double flux_gr=0;                ///< ground heat flux
        double L=0;                      ///< Obukhov length
        double zeta_m=0;                 ///< zm/L
        double zeta_s=0;                 ///< zs/L
        double zeta_o=0;                 ///< zo/L
        double zeta_t=0;                 ///< zt/L
        double sfc_T_new=0;              ///< surface temperature for next time step
        double sfc_q_new=0;              ///< surface moisture for next time step
        std::vector<double> soil_T_last; ///< soil temperature from previous time step
        std::vector<double> soil_q_last; ///< soil moisture from previous time step

        // Local atmospheric data
        double atm_U; ///< wind speed
        double atm_T; ///< temperature
        double atm_q; ///< mixing ratio
        double atm_p; ///< pressure
        double R_net; ///< net surface radiation
    
        // Local time data
        bool first=true;  ///< flag whether first time step or not
        int step_count=0; ///< number of times the LSM has been called
        double tstep=0;   ///< current time step
        double runtime=0; ///< current elapsed time
        
        // Local output information
        bool master;                     ///< flag whether UtahLSM column is first in a grid
        int j=0;                         ///< y-index of UtahLSM column if in a grid
        int i=0;                         ///< x-index of UtahLSM column if in a grid
        int output_counter=0;            ///< number of times output has been written
        std::vector<NcDim> dim_scalar_t; ///< dimension vector for time
        std::vector<NcDim> dim_scalar_z; ///< dimension vector for soil depth
        std::vector<NcDim> dim_vector;   ///< dimension vector for other fields
            
        /**
         * Struct to hold attributes of scalar fields.
         */
        struct AttScalar {
            double* data;                  ///< pointer to field
            std::string name;              ///< name of the field
            std::string long_name;         ///< description of the field
            std::string units;             ///< units of the field
            std::vector<NcDim> dimensions; ///< dimensions of the field
        };

        /**
         * Struct to hold attributes of vector fields.
         */
        struct AttVector {
            std::vector<double>* data;     ///< pointer to field
            std::string name;              ///< name of the field
            std::string long_name;         ///< description of the field
            std::string units;             ///< units of the field
            std::vector<NcDim> dimensions; ///< dimensions of the field
        };

        std::map<std::string,AttScalar> map_att_scalar; ///< map of scalar field name, attributes
        std::map<std::string,AttVector> map_att_vector; ///< map of vector field name, attributes
        std::vector<AttScalar> output_scalar;           ///< vector of scalar attributes to save
        std::vector<AttVector> output_vector;           ///< vector of vector attributes to save

        /**
         * Computes fluxes using similarity theory.
         *
         * @param[in] sfc_T surface temperature
         * @param[in] sfc_q surface mixing ratio
         */
        void computeFluxes(double sfc_T, double sfc_q);
        
        /**
         * Solves the surface energy budget.
         */
        void solveSEB();
        
        /**
         * Solves the surface moisture budget.
         */
        void solveSMB();
        
        /**
         * Solves the diffusion equation for soil heat.
         */
        void solveDiffusionHeat();

        /**
         * Solves the diffusion equation for soil moisture.
         */
        void solveDiffusionMois();
        
        /**
         * Solves the diffusion equation for soil heat and moisture.
         * 
         * @param[in] type 1=heat, 2=moisture
         */
        void solveDiffusion(int type);
        
        /**
         * Computes the surface energy budget.
         * 
         * @param[in] sfc_T surface temperature
         * @return surface energy budget
         */
        double computeSEB(double sfc_T);
        
        /**
         * Computes the derivative of the surface energy budget.
         * 
         * @param[in] sfc_T surface temperature
         * @return derivative of the surface energy budget
         */
        double computeDSEB(double sfc_T);
};

typedef void * LSMObject;      ///< Pointer representing UtahLSM object 
typedef void * SettingsObject; ///< Pointer representing Settings object 
typedef void * InputObject;    ///< Pointer representing Input object 
typedef void * OutputObject;   ///< Pointer representing Output object 

/** 
 * C-style interface for compatibility with other languages.
 */
extern "C" {
    /**
     * C-style wrapper for the UtahLSM constructor.
     *
     * @param[in]     settings Settings object
     * @param[in]     input Input object
     * @param[in]     output Output object
     * @param[in,out] ustar friction velocity
     * @param[in,out] flux_wT kinematic heat flux
     * @param[in,out] flux_wq kinematic moisture flux
     * @param[in]     j y-index of UtahLSM column if in a grid
     * @param[in]     i x-index of UtahLSM column if in a grid
     */
    LSMObject GetLSM(SettingsObject settings,InputObject input, OutputObject output,
                     double* ustar, double* flux_wT, 
                     double* flux_wq, int* j, int* i);
    
    /**
     * C-style wrapper for the updateFields function.
     *
     * @param[in] lsm UtahLSM object
     * @param[in] dt current time-step
     * @param[in] u wind speed
     * @param[in] T temperature
     * @param[in] q mixing ratio
     * @param[in] p pressure
     * @param[in] rad optional net radiation
     */
    void UpdateFields(LSMObject lsm,double* dt,double* u,double* T,double* q,double* p,double* rad);
    
    /**
     * C-style wrapper for the run function.
     * 
     * @param[in] lsm UtahLSM object
     */
    void Run(LSMObject lsm);

    /**
     * C-style wrapper for the save function.
     * 
     * @param[in] lsm UtahLSM object
     * @param[in] output Output object
     */
    void Save(LSMObject lsm, OutputObject output);
}

#endif