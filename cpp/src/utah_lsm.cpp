/*
 * UtahLSM
 * 
 * Copyright (c) 2017–2022 Jeremy A. Gibbs
 * Copyright (c) 2017–2022 Rob Stoll
 * Copyright (c) 2017–2022 Eric Pardyjak
 * Copyright (c) 2017–2022 Pete Willemsen
 * 
 * This file is part of UtahLSM.
 * 
 * This software is free and is distributed under the MIT License.
 * See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
 */

#include "utah_lsm.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>
#include <iomanip>
#include <span>

#include "constants.hpp"
#include "json.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "most.hpp"
#include "output.hpp"
#include "radiation.hpp"
#include "soil.hpp"

using json = nlohmann::json;

namespace {
    namespace c = constants;
}

// Constructor for UtahLSM class
UtahLSM :: UtahLSM(Input* input, Output* output, double& ustar, double& flux_wT,
                   double& flux_wq, int j, int i) : ustar(ustar),
                   flux_wT(flux_wT),flux_wq(flux_wq), j(j), i(i) {

    std::cout<<"[UtahLSM] \t Preparing to run"<<std::endl;

    // Input time section
    input->getItem(step_seb,"time","step_seb");
    input->getItem(step_dif,"time","step_dif");

    // Input grid section
    input->getItem(nx,"grid","nx");
    input->getItem(ny,"grid","ny");
    
    // Input length scale section
    input->getItem(z_o,"length","z_o");
    input->getItem(z_t,"length","z_t");
    input->getItem(z_m,"length","z_m");
    input->getItem(z_s,"length","z_s");
                 
    // Input soil section
    input->getItem(nsoilz,"soil","nsoil");
    input->getItem(soil_param,"soil","param");
    input->getItem(soil_model,"soil","model");
    input->getItem(soil_z,"soil","soil_z");
    input->getItem(soil_type,"soil","soil_type");
    input->getItem(soil_T,"soil","soil_T");
    input->getItem(soil_q,"soil","soil_q");
                       
    // Initialize new surface values for first run
    sfc_T_new = soil_T[0];
    sfc_q_new = soil_q[0];

    // Initialize history arrays for first run
    soil_T_last = soil_T;
    soil_q_last = soil_q;
    
    // Modify soil levels to be negative
    std::transform(soil_z.begin(), soil_z.end(), soil_z.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, -1.0));
    
    // Input radiation section
    input->getItem(comp_rad,"radiation","comp_rad");
    if (comp_rad==1) {
        input->getItem(utc_start,"radiation","utc_start");
        input->getItem(albedo,"radiation","albedo");
        input->getItem(emissivity,"radiation","emissivity");
        input->getItem(latitude,"radiation","latitude");
        input->getItem(longitude,"radiation","longitude");
        input->getItem(julian_day,"radiation","julian_day");
        
        // set initial time
        utc = utc_start;
        
        // convert latitude and longitude into radians
        latitude  = latitude * c::pi / 180.0;
        longitude = longitude * c::pi / 180.0;

        // Create radiation class
        radiation = new Radiation(latitude,longitude,albedo,emissivity);
    }
    // Create soil class
    soil = Soil::getModel(soil_type,soil_param,soil_model,nsoilz);

    // Input output section
    input->getItem(save_output, "output", "save");
    if (save_output) {
        
        // Get fields to save from user
        input->getItem(output_fields,"output","fields");
        if (output_fields[0]=="all") {
            output_fields.erase(output_fields.begin());
            output_fields = {"ust","shf","lhf","ghf","obl","soilt","soilq"};
        }
        
        // Decide if master UtahLSM column or not
        (j==0 && i==0) ? master=true : master=false;
        
        // Add dimensions
        NcDim t_dim, z_dim;
        if (master) {
            t_dim = output->addDimension("t");
            z_dim = output->addDimension("z",nsoilz);
        } else {
            t_dim = NcDim();
            z_dim = NcDim();
        }
        dim_scalar_t.push_back(t_dim);
        dim_scalar_z.push_back(z_dim);
        dim_vector.push_back(t_dim);
        dim_vector.push_back(z_dim);
        
        if (ny>1) {
            NcDim y_dim;
            if (master) {
                y_dim = output->addDimension("y",ny);
            }else {
                y_dim = NcDim();
            }
            dim_scalar_t.push_back(y_dim);
            dim_scalar_z.push_back(y_dim);
            dim_vector.push_back(y_dim);
        }
        if (nx>1) {
            NcDim x_dim;
            if (master) {
                x_dim = output->addDimension("x",nx);
            }else {
                x_dim = NcDim();
            }
            dim_scalar_t.push_back(x_dim);
            dim_scalar_z.push_back(x_dim);
            dim_vector.push_back(x_dim);
        }

        // Attributes for each field
        AttScalar att_time  = {&runtime, "time",  "time",               "s",     dim_scalar_t};
        AttScalar att_ust   = {&ustar,   "ust",   "friction velocity",  "m s-1", dim_scalar_t};
        AttScalar att_shf   = {&flux_sh, "shf",   "sensible heat flux", "W m-2", dim_scalar_t};
        AttScalar att_lhf   = {&flux_lh, "lhf",   "latent heat flux",   "W m-2", dim_scalar_t};
        AttScalar att_ghf   = {&flux_gr, "ghf",   "ground heat flux",   "W m-2", dim_scalar_t};
        AttScalar att_obl   = {&L,       "obl",   "Obukhov length",     "m",     dim_scalar_t};
        AttVector att_soilz = {&soil_z,  "soilz", "Obukhov length",     "m",     dim_scalar_z};
        AttVector att_soilt = {&soil_T,  "soilt", "soil temperature",   "K",     dim_vector};
        AttVector att_soilq = {&soil_q,  "soilq", "soil moisture",      "m3 m-3",dim_vector};
        
        // Map the name to attributes
        map_att_scalar.emplace("time", att_time);
        map_att_scalar.emplace("ust",  att_ust);
        map_att_scalar.emplace("shf",  att_shf);
        map_att_scalar.emplace("lhf",  att_lhf);
        map_att_scalar.emplace("ghf",  att_ghf);
        map_att_scalar.emplace("obl",  att_obl);
        map_att_vector.emplace("soilz",att_soilz);
        map_att_vector.emplace("soilt",att_soilt);
        map_att_vector.emplace("soilq",att_soilq);
        
        // We will always save time and depth
        output_scalar.push_back(map_att_scalar["time"]);
        output_vector.push_back(map_att_vector["soilz"]);
        
        // Create list of fields to save
        for (int i=0; i<output_fields.size(); i++) {
            std::string key = output_fields[i];
            if (map_att_scalar.count(key)) {
                output_scalar.push_back(map_att_scalar[key]);
            } else if (map_att_vector.count(key)) {
                output_vector.push_back(map_att_vector[key]);
            }
        }

        // Add scalar fields
        for (int i=0; i<output_scalar.size(); i++) {
            AttScalar att = output_scalar[i];
            if (master) {
                output->addField(att.name, att.units, att.long_name, att.dimensions);
            }
        }
        // Add vector fields
        for (int i=0; i<output_vector.size(); i++) {
            AttVector  att = output_vector[i];
            if (master) {
                output->addField(att.name, att.units, att.long_name, att.dimensions);
            }
        }
    }
}

// Update atmospheric quantities prior to solving
void UtahLSM :: updateFields(double dt,double u,double T,double q,double p,double rad=0) {
    
    tstep = dt;
    atm_U = u;
    atm_T = T;
    atm_q = q;
    atm_p = p;
    R_net = rad;
    runtime += tstep;
    
    // Run radiation model and update time/date if needed
    if (comp_rad==1) {
        utc = std::fmod(runtime,86400);
        julian_day += int(runtime/86400);
        R_net = radiation->computeNet(julian_day,utc,soil_T[0]);
    }
    
    // Keep winds from being exactly zero
    if (atm_U==0) atm_U = 0.1;
}

// Run UtahLSM
void UtahLSM :: run() {
    
    // Set initial new temp and moisture
    sfc_T_new = soil_T[0];
    sfc_q_new = soil_q[0];

    // Check if time to re-compute balances
    if ( (step_count % step_seb)==0 ) {
        solveSEB();
        solveSMB();
    } else {
        // just return new fluxes
        computeFluxes(soil_T[0],soil_q[0]);
    }

    // Save current temperature and moisture
    soil_T_last = soil_T;
    soil_q_last = soil_q;

    // check if time to compute diffusion
    if ( (step_count % step_dif)==0 ) {

        // Solve heat diffusion
        solveDiffusionHeat();

        // solve moisture diffusion
        solveDiffusionMois();
    }
    
    // Change flag of whether initial time
    if (first) first=false;

    // Increment step counter
    step_count += 1;
}

// Save current fields to the output file
void UtahLSM :: save(Output* output) {
    
    // Output size and location
    std::vector<size_t> scalar_index; 
    std::vector<size_t> scalar_size;
    std::vector<size_t> vector_index;
    std::vector<size_t> vector_size;
    std::vector<size_t> vector_index_z;
    std::vector<size_t> vector_size_z;
    if (ny==1 && nx==1) {
        scalar_index = {static_cast<unsigned long>(output_counter)};
        vector_index = {static_cast<size_t>(output_counter), 0};
        vector_size  = {1, static_cast<unsigned long>(nsoilz)};
    }
    else if (ny>1 && nx==1) {
        scalar_index   = {static_cast<unsigned long>(output_counter), static_cast<unsigned long>(j)};
        vector_index   = {static_cast<size_t>(output_counter), 0, static_cast<unsigned long>(j)};
        vector_size    = {1, static_cast<unsigned long>(nsoilz), 1};
        vector_index_z = {0, static_cast<unsigned long>(j)};
        vector_size_z  = {static_cast<unsigned long>(nsoilz),1};
    } else if (ny==1 && nx>1) {
        scalar_index   = {static_cast<unsigned long>(output_counter), static_cast<unsigned long>(i)};
        vector_index   = {static_cast<size_t>(output_counter), 0, static_cast<unsigned long>(i)};
        vector_size    = {1, static_cast<unsigned long>(nsoilz), 1};
        vector_index_z = {0, static_cast<unsigned long>(i)};
        vector_size_z  = {static_cast<unsigned long>(nsoilz),1};
    }
    else {
        scalar_index   = {static_cast<unsigned long>(output_counter), static_cast<unsigned long>(j),
                          static_cast<unsigned long>(i)};
        scalar_size    = {1, 1, 1};
        vector_index   = {static_cast<size_t>(output_counter), 0, static_cast<unsigned long>(j),
                          static_cast<unsigned long>(i)};
        vector_size    = {1, static_cast<unsigned long>(nsoilz),1, 1};
        vector_index_z = {0, static_cast<unsigned long>(j),static_cast<unsigned long>(i)};
        vector_size_z  = {static_cast<unsigned long>(nsoilz),1 ,1};
    }
    
    // Loop through scalar fields to save
    for (int i=0; i<output_scalar.size(); i++) {
        if (i==0 && master) {
            output->saveFieldScalar(output_scalar[i].name, scalar_index, output_scalar[i].data);
        } else if (i>0) {
            output->saveFieldScalar(output_scalar[i].name, scalar_index, output_scalar[i].data);
        }
    }
    // Loop through vector fields to save
    for (int i=0; i<output_vector.size(); i++) {
        
        // Soil depth is only saved once with no time component
        if (i==0 && output_counter==0) {
            if (ny==1 && nx==1)  {
                output->saveFieldVector(output_vector[i].name, *output_vector[i].data);
            } else {
                output->saveFieldVector(output_vector[i].name, vector_index_z,
                                        vector_size_z, *output_vector[i].data);
            }
        } else {
            output->saveFieldVector(output_vector[i].name, vector_index,
                                    vector_size, *output_vector[i].data);
        }
    }
    
    // Remove soil depth from fields to save after first time loop
    if (output_counter==0) output_vector.erase(output_vector.begin());
    
    // Increment for next time insertion
    output_counter +=1;
}

// Compute surface fluxes using Monin-Obukhov w/Dyer-Hicks
void UtahLSM :: computeFluxes(double sfc_T, double sfc_q) {
    
    // Local variables
    int max_iterations = 200;
    bool converged = false;
    double gnd_q, flux_wTv;
    double last_L, criteria = 0.1, ref_T = 300.;
    int int_depth;
    double heat_cap, dT, dz, K0, K1;
    
    // Compute surface mixing ratio
    gnd_q  = soil->surfaceMixingRatio(sfc_T,sfc_q,atm_p);
    // std::cout<<"computeFluxes, sfc mix ratio: "<<gnd_q<<std::endl;
    // std::cin.get();
    // Sensible flux, latent flux, ustar, and L
    for (int i=0; i<max_iterations; ++i) {
        
        // Compute ground flux
        // First time through we estimate based on Santanello and Friedl (2003)
        if ( (first)) {
            float A,B;
            if (soil_q[0]>=0.4) {
                A = 0.31;
                B = 74000;
            } else if (soil_q[0]<0.4 && soil_q[0] >= 0.25){
                A = 0.33;
                B = 85000;
            } else {
                A = 0.35;
                B = 100000;
            }
            flux_gr = R_net*A*std::cos((2*c::pi*(utc)+10800)/B);
        } else {
            double K0 = soil->conductivityThermal(soil_q[0],0);
            double K1 = soil->conductivityThermal(soil_q[1],1);
            double Kmid = 0.5*(K0 + K1);
            flux_gr = Kmid*(sfc_T - soil_T[1])/(soil_z[0]-soil_z[1]);
        }
        
        // Compute friction velocity
        ustar = atm_U*most::fm(z_m/z_o,zeta_m,zeta_o);
        
        // Compute heat flux
        flux_wT = (sfc_T-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // Compute latent flux
        if ( (first) && (i == 0)) {
            flux_wq = (R_net - flux_gr - flux_wT*c::rho_air*c::Cp_air)/(c::rho_air*c::Lv);
            gnd_q = atm_q + flux_wq / (ustar*most::fh(z_s/z_t,zeta_s,zeta_t));
            soil_q[0] = soil->surfaceWaterContentEstimate(soil_T[0],gnd_q, atm_p);
            sfc_q_new = soil_q[0];
        } else {
            flux_wq = (gnd_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        }
        
        // Compute virtual heat flux
        flux_wTv = flux_wT + ref_T*0.61*flux_wq;
        
        // Compute L
        last_L = L;
        L      = -std::pow(ustar,3.)*ref_T/(c::vonk*c::grav*flux_wTv);
        
        // Bounds check on L
        if (z_m/L > 5.)  L = z_m/5.;
        if (z_m/L < -5.) L = -z_m/5.;
        
        // Update zeta terms
        zeta_m = z_m/L;
        zeta_s = z_s/L;
        zeta_o = z_o/L;
        zeta_t = z_t/L;
        
        // Check for convergence
        converged = std::abs(last_L-L) <= criteria;
        if (converged) {
            flux_sh = c::rho_air*c::Cp_air*flux_wT;
            flux_lh = c::rho_air*c::Lv*flux_wq;
            break;
        }
    }
    
    if (!converged) {
        std::cout<<std::endl;
        std::cout<<"[Fluxes] \t Converge failed"<<std::endl;
        std::cout<<L<<" "<<flux_wq<<std::endl;
        throw(1);
    }
}

// Solve the surface energy balance
void UtahLSM :: solveSEB() {
    
    // Local variables
    int max_iter_temp = 200;
    int max_iter_flux = 200;
    double dTs, dTs_old;
    double SEB, dSEB_dT, SEB_l, SEB_h;
    double temp_l, temp_h, last_T, last_F;
    double temp_1 = soil_T[0] - 1;
    double temp_2 = soil_T[0] + 1;
    double temp_criteria = 0.001;
    double flux_criteria = 0.001;
    
    // Compute SEB using current bracketed temperatures
    SEB_l = computeSEB(temp_1);
    SEB_h = computeSEB(temp_2);

    // Dynamic bracket adjustments
    bool out_of_bracket = (SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0);
    
    while (out_of_bracket) {
        
        // Expand brackets by 1 K
        temp_1 -= 1;
        temp_2 += 1;
        
        // Recompute SEB at brackets
        SEB_l = computeSEB(temp_1);
        SEB_h = computeSEB(temp_2);
        
        // Check for proper brackets
        out_of_bracket = (SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0);
    }
    
    if ((SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0)) {
        throw("Please adjust brackets for Ts");
    }
    
    // If SEB from low bracket Ts = 0, then that value of Ts is solution
    if (SEB_l == 0.0) sfc_T_new = temp_1;
    
    // If SEB from high bracket Ts = 0, then that value of Ts is solution
    if (SEB_h == 0.0) sfc_T_new = temp_2;
    
    // Orient the solutions such that SEB(temp_l) < 0;
    if (SEB_l < 0.0) {
        temp_l = temp_1;
        temp_h = temp_2;
    } else {
        temp_l = temp_2;
        temp_h = temp_1;
    }
    
    // Prepare for convergence looping
    dTs = std::abs(temp_h-temp_l);
    
    // Convergence loop for flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // Convergence loop for temperature
        for (int tt = 0; tt < max_iter_temp; tt++) {
            
            // Compute SEB and dSEB_dTs
            SEB     = computeSEB(sfc_T_new);
            dSEB_dT = computeDSEB(sfc_T_new);
            
            // Update brackets
            if (SEB<0.) temp_l = sfc_T_new;
            if (SEB>0.) temp_h = sfc_T_new;
            
            // Bracket and bisect temperature if Newton out of range
            if ((((sfc_T_new-temp_h)*dSEB_dT-SEB)*((sfc_T_new-temp_l)*dSEB_dT-SEB)>0.0)
                || (std::abs(2*SEB) > std::abs(dTs_old*dSEB_dT))) {
                dTs_old   = dTs;
                dTs       = 0.5*(temp_h-temp_l);
                last_T    = sfc_T_new;
                sfc_T_new = temp_l + dTs;
                if (temp_l == sfc_T_new) break;
            } else {
                dTs_old   = dTs;
                dTs       = SEB / dSEB_dT;
                last_T    = sfc_T_new;
                sfc_T_new = sfc_T_new - dTs;
                if (last_T == sfc_T_new) break;
            }
            
            // Check for convergence
            if (std::abs( (sfc_T_new-last_T)/last_T) <= temp_criteria) break;
            
            // If convergence fails, recompute flux
            // computeFluxes(soil_T[0],soil_q[0]);
        }
        
        // Save current flux for convergence criteria
        last_F = flux_wT;
        
        // Recompute heat flux using new temperature
        computeFluxes(sfc_T_new,sfc_q_new);
        
        // Check for convergence
        if (std::abs(flux_wT-last_F) <= flux_criteria) {
            double Qh = c::rho_air*c::Cp_air*flux_wT;
            double Ql = c::rho_air*c::Lv*flux_wq;
            double Qg = flux_gr;
            break;
        }
        
        // If flux fails to converge, split temperature difference
        sfc_T_new = 0.5*(sfc_T_new + last_T);
    }
}

// Compute surface energy budget
double UtahLSM :: computeSEB(double sfc_T) {
    
    // Local variables
    double Qg, Qh, Ql, SEB;
    
    // Compute fluxes using passed in values
    computeFluxes(sfc_T,sfc_q_new);
    
    // Write sensible and latent heat fluxes in [W/m^2]
    Qh = c::rho_air*c::Cp_air*flux_wT;
    Ql = c::rho_air*c::Lv*flux_wq;
    Qg = flux_gr;
    
    // Compute surface energy balance
    SEB = R_net - Qg - Qh - Ql;
    return SEB;
}

// Computes derivative of surface energy budget wrt surface temperature
double UtahLSM :: computeDSEB(double sfc_T) {
    
    // Local variables
    double heat_cap, dSEB_dT;
    
    // Compute derivative of SEB wrt temperature
    heat_cap = soil->heatCapacity(sfc_q_new,0);
    dSEB_dT = 4.*emissivity*c::sb*std::pow(sfc_T,3.)
    + c::rho_air*c::Cp_air*ustar*most::fh(z_s/z_t,zeta_s,zeta_t)
    + heat_cap*(soil_z[0]-soil_z[1])/(2*tstep);
    return dSEB_dT;
}

// Solve the surface moisture balance
void UtahLSM :: solveSMB() {
    
    // Local variables
    bool converged;
    int max_iter_flux = 200;
    double E,flux_sm_last, flux_sm, flux_sm2;;
    double psi0, psi1, K0, K1, K_avg, D0, D1, D_avg;
    double delta = 0.5, flux_criteria = .001; 
    
    // Moisture potential at first two levels below ground
    psi0 = soil->waterPotential(soil_q[0], 0);
    psi1 = soil->waterPotential(soil_q[1], 1);
    
    // Compute initial soil moisture flux
    K0    = soil->conductivityMoisture(soil_q[0],0);
    K1    = soil->conductivityMoisture(soil_q[1],1);
    K_avg = 0.5*(K0+K1);
    D0    = soil->diffusivityMoisture(soil_q[0],0);
    D1    = soil->diffusivityMoisture(soil_q[1],1);
    D_avg = 0.5*(D0+D1);
    flux_sm  = c::rho_wat*K_avg*((psi0 - psi1)/(soil_z[0]-soil_z[1]) + 1);
    flux_sm = c::rho_wat*D_avg*(soil_q[0]-soil_q[1])/(soil_z[0]-soil_z[1])
               + c::rho_wat*K_avg;
    
    // Compute evaporation
    E = c::rho_air*flux_wq;
    
    // Convergence loop for moisture flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // Save soil moisture flux for convergence test
        flux_sm_last = flux_sm;
        
        // Compute new weighted soil moisture flux
        flux_sm = delta*flux_sm_last - (1.-delta)*E;
        
        // Re-compute moisture potential
        psi0    = psi1 + (soil_z[0]-soil_z[1])*((flux_sm/(c::rho_wat*K_avg))-1);
        if (psi0 > soil->properties[0]->psi_sat) {
            psi0 = soil->properties[0]->psi_sat;
        }
        
        // Update soil moisture
        sfc_q_new = soil->surfaceWaterContent(psi0);
        double gnd_q  = soil->surfaceMixingRatio(sfc_T_new,sfc_q_new,atm_p);
        E = c::rho_air*(gnd_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // Update soil moisture transfer
        K0    = soil->conductivityMoisture(sfc_q_new,0);
        K1    = soil->conductivityMoisture(sfc_q_new,1);
        K_avg = 0.5*(K0+K1);
        
        // Check for convergence
        converged = std::abs((E + flux_sm)/E) <=flux_criteria;

        if (converged) {
            break;
        }
    }
}

void UtahLSM :: solveDiffusionHeat() {
    
    // Local variables
    double AB  = 1.0;
    double AF  = 1.0-AB;
    double dz  = soil_z[0] - soil_z[1];
    double dz2 = std::pow(dz,2);
    double Cp, Cm, CBp, CBm, CB, CFp, CFm, CF;

    std::vector<double> K(nsoilz,0.0);
    std::vector<double> K_mid(nsoilz-1,0.0);
    std::vector<double> z_mid(nsoilz-1,0.0);
    std::vector<double> r(nsoilz-1,0.0);
    std::vector<double> e(nsoilz-1,0.0);
    std::vector<double> f(nsoilz-1,0.0);
    std::vector<double> g(nsoilz-1,0.0);
    std::cout<<"IN solveDiffusionHeat"<<std::endl;
    for (int i=0; i<nsoilz-1; i++) {

        K[i]     = soil->diffusivityThermal(soil_q[i],i);
        K[i+1]   = soil->diffusivityThermal(soil_q[i+1],i+1);
        K_mid[i] = 0.5*(K[i]+K[i+1]);
        z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
    }

    // Get the time step restriction
    double Kmax, dt_T;
    dt_T = 1;

    // Loop through diffusion by substep
    for (int t=0; t<=tstep; t+=dt_T) {

        // Set up and solve a tridiagonal matrix
        // AT(n+1) = r(n), where n denotes the time level
        // e, f, g the components of A matrix
        // T(n+1)  the soil temperature vector at t=n+1
        // r(n)    the soil temperature vector at t=n multiplied by coefficients

        // Matrix coefficients for first level below surface
        Cp  = step_dif * dt_T * K_mid[0] / dz2;
        Cm  = step_dif * dt_T * K_mid[1] / dz2;
        CBp = -AB * Cp;
        CBm = -AB * Cm;
        CB  = 1 - CBp - CBm;
        CFp = AF * Cp;
        CFm = AF * Cm;
        CF  = 1 - CFp - CFm;

        e[0] = 0;
        f[0] = CB;
        g[0] = CBm;
        r[0] = CFp * soil_T[0] + CF * soil_T[1] + CFm * soil_T[2] - CBp * sfc_T_new;
        
        // Matrix coefficients for the interior levels
        for (int i=1; i<nsoilz-2; i++) {

            // for soil_T in this loop:
            // i   -> j+1 level
            // i+1 -> j   level
            // i+2 -> j-1 level
            Cp  = step_dif * dt_T * K_mid[i] / dz2;
            Cm  = step_dif * dt_T * K_mid[i+1] / dz2;
            CBp = -AB * Cp;
            CBm = -AB * Cm;
            CB  = 1 - CBp - CBm;
            CFp = AF * Cp;
            CFm = AF * Cm;
            CF  = 1 - CFp - CFm;

            e[i] = CBp;
            f[i] = CB;
            g[i] = CBm;
            r[i] = CFp * soil_T[i] + CF * soil_T[i+1] + CFm * soil_T[i+2];
        }

        // Matrix coefficients for bottom level
        int j = nsoilz-2;

        Cp  = step_dif * dt_T * K_mid[j] / dz2;
        Cm  = step_dif * dt_T * K_mid[j] / dz2;
        CBp = -AB * Cp;
        CBm = -AB * Cm;
        CB  = 1 - CBp - CBm;
        CFp = AF * Cp;
        CFm = AF * Cm;
        CF  = 1 - CFp - CFm;

        e[j] = (CBp - CBm);
        f[j] = (CB + 2 * CBm);
        g[j] = 0;
        r[j] = (CFp - CFm) * soil_T[j] + (CF + 2* CFm) * soil_T[j+1];

        // now we can add new sfc T to column array
        soil_T[0] = sfc_T_new;

        // Solve the tridiagonal system
        try {
            // we only need to send the layers below surface
            std::span<double> subsfc_T(soil_T.cbegin() + 1, soil_T.size() - 1);
            matrix::tridiagonal(e,f,g,r,subsfc_T);
        } catch(std::string &e) {
            std::cout<<e<<std::endl;
            std::exit(0);
        }

        // solve K for this step to get a new dt
        for (int i=0; i<nsoilz-1; i++) {
            K[i]     = soil->diffusivityThermal(soil_q[i],i);
            K[i+1]   = soil->diffusivityThermal(soil_q[i+1],i+1);
            K_mid[i] = 0.5*(K[i]+K[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }

        // Get the time step restriction
        Kmax = *std::max_element(K.begin(), K.end());
        dt_T = dz2 / (2*Kmax); 
    }
}

void UtahLSM :: solveDiffusionMois() {
    
    // Local variables
    double AB  = 1.0;
    double AF  = 1.0-AB;
    double dz  = soil_z[0] - soil_z[1];
    double dz2 = std::pow(dz,2);
    double Cp,Cm;
    double Cpd, Cmd, Cpk, Cmk;
    double CBpd,CBmd,CBpk,CBmk;
    double CFpd,CFmd,CFpk,CFmk;
    double CBp, CBm, CB, CFp, CFm, CF;

    std::vector<double> K_lin(nsoilz,0.0);
    std::vector<double> D(nsoilz,0.0);
    std::vector<double> D_mid(nsoilz-1,0.0);
    std::vector<double> z_mid(nsoilz-1,0.0);
    std::vector<double> r(nsoilz-1,0.0);
    std::vector<double> e(nsoilz-1,0.0);
    std::vector<double> f(nsoilz-1,0.0);
    std::vector<double> g(nsoilz-1,0.0);

    // Get the time step restriction
    double Dmax, dt_q;
    dt_q = 1;
    
    // Loop through diffusion by substep
    for (int t=0; t<=tstep; t+=dt_q) {

        // Set up and solve a tridiagonal matrix
        // AT(n+1) = r(n), where n denotes the time level
        // e, f, g the components of A matrix
        // T(n+1)  the soil temperature vector at t=n+1
        // r(n)    the soil temperature vector at t=n multiplied by coefficients

        // first soil level below the surface
        // common coefficients
        Cpd  = step_dif * dt_q * D_mid[0] / dz2;
        Cmd  = step_dif * dt_q * D_mid[1] / dz2;
        Cpk  = step_dif * dt_q * K_lin[0] / (2*dz);
        Cmk  = step_dif * dt_q * K_lin[2] / (2*dz);

        // coefficients for backward scheme
        CBpd = -AB * Cpd;
        CBmd = -AB * Cmd;
        CBpk = -AB * Cpk;
        CBmk = -AB * Cmk;
        CB   = (1 - CBpd - CBmd);
        CBp  = CBpd + CBpk;
        CBm  = CBmd - CBmk;

        // coefficients for forward scheme
        CFpd = AF * Cpd;
        CFmd = AF * Cmd;
        CFpk = AF * Cpk;
        CFmk = AF * Cmk;
        CF   = (1 - CFpd - CFmd);
        CFp  = CFpd + CFpk;
        CFm  = CFmd - CFmk;

        // matrix components
        e[0] = 0;
        f[0] = CB;
        g[0] = CBm;
        r[0] = CFp * soil_q[0] + CF * soil_q[1] + CFm * soil_q[2] - CBp * sfc_q_new;

        // interior soil levels
        for (int i=1; i<nsoilz-2; i++) {

            // for soil_T in this loop:
            // i   -> j+1 level
            // i+1 -> j   level
            // i+2 -> j-1 level

            // common coefficients
            Cpd  = step_dif * dt_q * D_mid[i] / dz2;
            Cmd  = step_dif * dt_q * D_mid[i+1] / dz2;
            Cpk  = step_dif * dt_q * K_lin[i] / (2*dz);
            Cmk  = step_dif * dt_q * K_lin[i+2] / (2*dz);

            // coefficients for backward scheme
            CBpd = -AB * Cpd;
            CBmd = -AB * Cmd;
            CBpk = -AB * Cpk;
            CBmk = -AB * Cmk;
            CB   = (1 - CBpd - CBmd);
            CBp  = CBpd + CBpk;
            CBm  = CBmd - CBmk;

            // coefficients for forward scheme
            CFpd = AF * Cpd;
            CFmd = AF * Cmd;
            CFpk = AF * Cpk;
            CFmk = AF * Cmk;
            CF   = (1 - CFpd - CFmd);
            CFp  = CFpd + CFpk;
            CFm  = CFmd - CFmk;

            // matrix components
            e[i] = CBp;
            f[i] = CB;
            g[i] = CBm;
            r[i] = CFp * soil_q[i] + CF * soil_q[i+1] + CFm * soil_q[i+2];
        }

        // Matrix coefficients for bottom level
        int j = nsoilz-2;

        // common coefficients
        Cpd  = step_dif * dt_q * D_mid[j] / dz2;
        Cmd  = step_dif * dt_q * D_mid[j] / dz2;
        Cpk  = step_dif * dt_q * K_lin[j] / (2*dz);
        Cmk  = step_dif * dt_q * K_lin[j] / (2*dz);

        // coefficients for backward scheme
        CBpd = -AB * Cpd;
        CBmd = -AB * Cmd;
        CBpk = -AB * Cpk;
        CBmk = -AB * Cmk;
        CB   = (1 - CBpd - CBmd);
        CBp  = CBpd + CBpk;
        CBm  = CBmd - CBmk;

        // coefficients for forward scheme
        CFpd = AF * Cpd;
        CFmd = AF * Cmd;
        CFpk = AF * Cpk;
        CFmk = AF * Cmk;
        CF   = (1 - CFpd - CFmd);
        CFp  = CFpd + CFpk;
        CFm  = CFmd - CFmk;

        // matrix components
        e[j] = (CBp - CBm);
        f[j] = (CB + 2 * CBm);
        g[j] = 0;
        r[j] = (CFp - CFm) * soil_q[j] + (CF + 2 * CFm) * soil_q[j+1];

        // now we can add new sfc q to column array
        soil_q[0] = sfc_q_new;

        // Solve the tridiagonal system
        try {
            // we only need the layers below the surface
            std::span<double> subsfc_q(soil_q.data() + 1, soil_q.size() - 1);
            matrix::tridiagonal(e,f,g,r,subsfc_q);
        } catch(std::string &e) {
            std::cout<<e<<std::endl;
            std::exit(0);
        }

        // solve for D to get new dt
        for (int i=0; i<nsoilz-1; i++) {
            D[i]     = soil->diffusivityMoisture(soil_q[i],i);
            D[i+1]   = soil->diffusivityMoisture(soil_q[i+1],i+1);
            D_mid[i] = 0.5*(D[i]+D[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);

            // linearized K
            K_lin[i] = soil->conductivityMoisture(soil_q[i],i)/soil_q[i];
            if (i==nsoilz-2) {
                K_lin[i+1] = soil->conductivityMoisture(soil_q[i+1],i+1)/soil_q[i+1];
            }
        }

        // get new dt
        Dmax = *std::max_element(D.begin(), D.end());
        dt_q = dz2 / (2*Dmax); 
    }
}

//////////////////////////////////////////////////////////////
// C-style interface for compatibility with other languages //
//////////////////////////////////////////////////////////////

// C-style wrapper for the UtahLSM constructor
LSMObject GetLSM(InputObject input, OutputObject output,
                 double* ustar, double* flux_wT, 
                 double* flux_wq, int* j, int* i) {
    
    // Get input and output objects
    Input* input_obj = (Input*)input;
    Output* output_obj = (Output*)output;
    
    // Create lsm object
    UtahLSM* lsm = new UtahLSM(input_obj,output_obj,*ustar,*flux_wT,*flux_wq,*j,*i);

    return (LSMObject)lsm;
}

// C-style wrapper for the updateFields function
void UpdateFields(LSMObject lsm,double* dt,double* u,double* T,double* q,double* p,double* rad=0) {
        
    // Get LSM object
    UtahLSM* lsm_obj = (UtahLSM*)lsm;
    
    // Update fields
    lsm_obj->updateFields(*dt,*u,*T,*q,*p,*rad);
    
    return;
    
}

// C-style wrapper for the run function
void Run(LSMObject lsm) {
    
    // Get LSM object
    UtahLSM* lsm_obj = (UtahLSM*)lsm;
    
    // Run the lsm
    lsm_obj->run();
}

// C-style wrapper for the save function
void Save(LSMObject lsm, OutputObject output) {
        
    // Get LSM and output objects
    UtahLSM* lsm_obj = (UtahLSM*)lsm;
    Output* output_obj = (Output*)output;
    
    // Save data to file
    lsm_obj->save(output_obj);
}