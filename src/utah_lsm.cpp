/*
 * UtahLSM
 * 
 * Copyright (c) 2019 Jeremy A. Gibbs
 * Copyright (c) 2019 Pete Willemsen
 * Copyright (c) 2019 Rob Stoll
 * Copyright (c) 2019 Eric Pardyjak
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
        AttScalar att_shf   = {&flux_wT, "shf",   "sensible heat flux", "W m-2", dim_scalar_t};
        AttScalar att_lhf   = {&flux_wq, "lhf",   "latent heat flux",   "W m-2", dim_scalar_t};
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
    
    // Save previoius temp and moisture
    surf_T_last = soil_T[0];
    surf_q_last = soil_q[0];
    
    // Check if time to re-compute balances
    if ( (step_count % step_seb)==0 ) {

        // Solve the surface energy balance
        solveSEB();

        // Solve the surface moisture balance
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
        solveDiffusion(1);

        // solve moisture diffusion
        solveDiffusion(2);
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
    int depth = nsoilz;
    int int_depth;
    double heat_cap, dT, dz, K0, K1;
    
    // Compute surface mixing ratio
    gnd_q  = soil->surfaceMixingRatio(sfc_T,sfc_q,atm_p);
    
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
            // Compute heat flux within soil and find depth of minimum
            // This is the depth to integrate time change of T
            double hfs;
            double min_sflux = 10000.;
            for (int d=0; d<nsoilz-1; ++d) {
                K0 = soil->conductivityThermal(soil_q[d],d);
                K1 = soil->conductivityThermal(soil_q[d+1],d+1);
                hfs = (d==0) ? 0.5*(K0+K1) * (sfc_T - soil_T[d+1])/(soil_z[d]-soil_z[d+1]):
                0.5*(K0+K1) * (soil_T[d] - soil_T[d+1])/(soil_z[d]-soil_z[d+1]);
                if (std::abs(hfs)<std::abs(min_sflux)) {
                    min_sflux = hfs;
                    int_depth = d+1;
                }
            }
            // Integrate time change to depth of minimum flux
            flux_gr = 0;
            for (int d=0; d<int_depth; ++d) {
                
                heat_cap = soil->heatCapacity(soil_q[d],d);
                dT = (d==0) ? sfc_T-soil_T_last[d] : soil_T[d]-soil_T_last[d];
                if (d==0) {
                    dz = soil_z[d] - soil_z[d+1];
                } else if (d==int_depth-1){
                    dz = soil_z[d-1] - soil_z[d];
                } else {
                    dz = soil_z[d-1] - soil_z[d+1];
                }
                flux_gr = flux_gr + 0.5*( heat_cap*dT*dz/tstep);
            }
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
            break;
        }
    }
    
    if (!converged) {
        std::cout<<std::endl;
        std::cout<<"[Fluxes] \t We did not converge"<<std::endl;
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
        
        // Expand brackets by 5 K
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
    if (SEB_l == 0.0) soil_T[0] = temp_1;
    
    // If SEB from high bracket Ts = 0, then that value of Ts is solution
    if (SEB_h == 0.0) soil_T[0] = temp_2;
    
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
            SEB     = computeSEB(soil_T[0]);
            dSEB_dT = computeDSEB(soil_T[0]);
            
            // Update brackets
            if (SEB<0.) temp_l = soil_T[0];
            if (SEB>0.) temp_h = soil_T[0];
            
            // Bracket and bisect temperature if Newton out of range
            if ((((soil_T[0]-temp_h)*dSEB_dT-SEB)*((soil_T[0]-temp_l)*dSEB_dT-SEB)>0.0)
                || (std::abs(2*SEB) > std::abs(dTs_old*dSEB_dT))) {
                dTs_old   = dTs;
                dTs       = 0.5*(temp_h-temp_l);
                last_T    = soil_T[0];
                soil_T[0] = temp_l + dTs;
                if (temp_l == soil_T[0]) break;
            } else {
                dTs_old   = dTs;
                dTs       = SEB / dSEB_dT;
                last_T    = soil_T[0];
                soil_T[0] = soil_T[0] - dTs;
                if (last_T == soil_T[0]) break;
            }
            
            // Check for convergence
            if (std::abs( (soil_T[0]-last_T)/last_T) <= temp_criteria) break;
            
            // If convergence fails, recompute flux
            computeFluxes(soil_T[0],soil_q[0]);
        }
        
        // Save current flux for convergence criteria
        last_F = flux_wT;
        
        // Recompute heat flux using new temperature
        computeFluxes(soil_T[0],soil_q[0]);
        
        // Check for convergence
        if (std::abs(flux_wT-last_F) <= flux_criteria) {
            break;
        }
        
        // If flux fails to converge, split temperature
        soil_T[0] = 0.5*(soil_T[0] + last_T);
        
        // Recompute stability corrections and fluxes
        computeFluxes(soil_T[0],soil_q[0]);
    }
}

// Compute surface energy budget
double UtahLSM :: computeSEB(double sfc_T) {
    
    // Local variables
    double Qg, Qh, Ql, SEB, K;
    
    // Compute fluxes using passed in values
    computeFluxes(sfc_T,soil_q[0]);
    
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
    heat_cap = soil->heatCapacity(soil_q[0],0);
    dSEB_dT = - 4.*emissivity*c::sb*std::pow(sfc_T,3.)
    - c::rho_air*c::Cp_air*ustar*most::fh(z_s/z_t,zeta_s,zeta_t)
    - heat_cap*(soil_z[0]-soil_z[1])/(2*tstep);
    return dSEB_dT;
}

// Solve the surface energy balance
void UtahLSM :: solveSMB() {
    
    // Local variables
    bool converged;
    int max_iter_flux = 200;
    double E,flux_sm_last, flux_sm, flux_sm2;;
    double psi0, psi1, K0, K1, K_avg, D0, D1, D_avg;
    double delta = 0.8, flux_criteria = .001; 
    
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
    flux_sm2 = c::rho_wat*D_avg*(soil_q[0]-soil_q[1])/(soil_z[0]-soil_z[1])
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
        soil_q[0] = soil->surfaceWaterContent(psi0);

        double gnd_q  = soil->surfaceMixingRatio(soil_T[0],soil_q[0],atm_p);

        E = c::rho_air*(gnd_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // Update soil moisture transfer
        K0    = soil->conductivityMoisture(soil_q[0],0);
        K1    = soil->conductivityMoisture(soil_q[1],1);
        K_avg = 0.5*(K0+K1);
        
        // Check for convergence
        converged = std::abs((E + flux_sm)/E) <=flux_criteria;
    }
}

// Integrate soil heat diffusion equation
void UtahLSM :: solveDiffusion(int type) {
    
    // Local variables
    double surf_scalar_last;
    double dKdz, C_tp, C_tm, C_p, C_m;
    double dz_p, dz_m, a_m, a_p, a_o;
    std::vector<double> K(nsoilz,0.0);
    std::vector<double> D(nsoilz,0.0);
    std::vector<double> K_mid(nsoilz-1,0.0);
    std::vector<double> z_mid(nsoilz-1,0.0);
    std::vector<double> r(nsoilz-1,0.0);
    std::vector<double> e(nsoilz-1,0.0);
    std::vector<double> f(nsoilz-1,0.0);
    std::vector<double> g(nsoilz-1,0.0);
    std::vector<double> scalar(nsoilz);
    
    // Set up and solve a tridiagonal matrix
    // AT(n+1) = r(n), where n denotes the time level
    // e, f, g the components of A matrix
    // T(n+1)  the soil temperature vector at t=n+1
    // r(n)    the soil temperature vector at t=n multiplied by coefficients
    
    // Interpolate soil_z and D_n to mid-points
    if (type==1) {
        for (int i=0; i<nsoilz-1; i++) {
            K[i]     = soil->conductivityThermal(soil_q[i],i);
            K[i+1]   = soil->conductivityThermal(soil_q[i+1],i);
            D[i]     = soil->diffusivityThermal(K[i],soil_q[i],i);
            D[i+1]   = soil->diffusivityThermal(K[i+1],soil_q[i+1],i);
            K_mid[i] = 0.5*(D[i]+D[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar = soil_T;
        surf_scalar_last = soil_T_last[0];
    } else {
        for (int i=0; i<nsoilz-1; i++) {
            K[i]       = soil->conductivityMoisture(soil_q[i],i);
            K[i+1]     = soil->conductivityMoisture(soil_q[i+1],i);
            D[i]       = soil->diffusivityMoisture(soil_q[i],i);
            D[i+1]     = soil->diffusivityMoisture(soil_q[i+1],i+1);
            K_mid[i]   = 0.5*(D[i]+D[i+1]);
            z_mid[i]   = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar           = soil_q;
        surf_scalar_last = soil_q_last[0];
    }
    
    // Matrix coefficients for first level below surface
    C_p  = tstep*K_mid[0]/(2*(z_mid[0]-z_mid[1])*(soil_z[0]-soil_z[1]));
    C_m  = tstep*K_mid[1]/(2*(z_mid[0]-z_mid[1])*(soil_z[1]-soil_z[2]));
    C_tp = (1.0 + C_p + C_m);
    C_tm = (1.0 - C_p - C_m);
    
    f[0] =  C_tp;
    g[0] = -C_m;
    r[0] =  C_p*surf_scalar_last + C_tm*scalar[1] + C_m*scalar[2] + C_p*scalar[0];
    
    // For moisture, we need additional dKn/dz term
    // use 3-pt stencil that allows non-unform spacing
    if (type==2) {
        
        dz_p = soil_z[0] - soil_z[1];
        dz_m = soil_z[1] - soil_z[2];
        a_m  = -dz_p / (dz_m * (dz_p + dz_m));
        a_o  = (dz_p - dz_m) / (dz_p * dz_m);
        a_p  = dz_m / (dz_p * (dz_p + dz_m) );
        
        dKdz = tstep*(K[0]*a_p + K[1]*a_o + K[2]*a_m);
        r[0] = r[0] + dKdz;
    }
    
    // Matrix coefficients for the interior levels
    for (int i=1; i<nsoilz-2; i++) {
        C_p  = tstep*K_mid[i]  / (2*(z_mid[i]-z_mid[i+1])*(soil_z[i]  -soil_z[i+1]));
        C_m  = tstep*K_mid[i+1]/ (2*(z_mid[i]-z_mid[i+1])*(soil_z[i+1]-soil_z[i+2]));
        C_tp = (1 + C_p + C_m);
        C_tm = (1 - C_p - C_m);
        
        e[i] = -C_p;
        f[i] =  C_tp;
        g[i] = -C_m;
        r[i] =  C_p*scalar[i] + C_tm*scalar[i+1] + C_m*scalar[i+2];
        
        // For moisture, we need additional dKn/dz term
        // use 3-pt stencil that allows non-unform spacing
        if (type==2) {
            dz_p = soil_z[i] - soil_z[i+1];
            dz_m = soil_z[i+1] - soil_z[i+2];
            a_m  = -dz_p / (dz_m * (dz_p + dz_m));
            a_o  = (dz_p - dz_m) / (dz_p * dz_m);
            a_p  = dz_m / (dz_p * (dz_p + dz_m) );
            dKdz = tstep*(K[i]*a_p + K[i+1]*a_o + K[i+2]*a_m);
            r[i] = r[i] + dKdz;
        }
    }
    
    // Matrix coefficients for bottom level
    // first, construct ghost values
    int j       = nsoilz-2;
    double z_g  = 2*soil_z[j+1] - soil_z[j];
    double z_mg = (soil_z[j+1] + z_g) / 2.;
    
    C_p  = tstep*K_mid[j]/(2*(z_mid[j]-z_mg)*(soil_z[j]-soil_z[j+1]));
    C_m  = tstep*K_mid[j]/(2*(z_mid[j]-z_mg)*(soil_z[j+1]-z_g));
    C_tp = (1 + C_p + C_m);
    C_tm = (1 - C_p - C_m);
    
    e[j] = C_m  - C_p;
    f[j] = C_tp - 2* C_m;
    r[j] = (C_p - C_m)*scalar[j] + (C_tm+2*C_m)*scalar[j+1];
    
    // For moisture, we need additional dKn/dz term
    // use simple 2-pt backward Euler approximation
    if (type==2) {
        dKdz = tstep*(K[nsoilz-2]-K[nsoilz-1])/(soil_z[nsoilz-2]-soil_z[nsoilz-1]);
        r[nsoilz-2] = r[nsoilz-2] + dKdz;
    }
    
    // Solve the tridiagonal system
    try {
        if (type==1) matrix::tridiagonal(e,f,g,r,soil_T);
        if (type==2) matrix::tridiagonal(e,f,g,r,soil_q);
    } catch(std::string &e) {
        std::cout<<e<<std::endl;
        std::exit(0);
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