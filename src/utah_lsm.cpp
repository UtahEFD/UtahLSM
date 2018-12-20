//
//  utah_lsm.cpp
//  
//  This class handles the UtahLSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include "utah_lsm.hpp"
#include "json.hpp"
#include "constants.hpp"
#include "input.hpp"
#include "output.hpp"
#include "soil.hpp"
#include "radiation.hpp"
#include "most.hpp"
#include "matrix.hpp"
#include<fstream>

using json = nlohmann::json;

namespace {
    namespace c = Constants;
}

///////////////////////
// Private functions //
///////////////////////

UtahLSM :: UtahLSM(Input* input, Output* output, double& ustar, double& flux_wT,
                   double& flux_wq, int j, int i) : ustar(ustar),
                   flux_wT(flux_wT),flux_wq(flux_wq), j(j), i(i) {

    std::cout<<"[UtahLSM] \t Preparing to run"<<std::endl;
    
    // grid section
    input->getItem(nx,"grid","nx");
    input->getItem(ny,"grid","ny");
    
    // length scale section
    input->getItem(z_o,"length","z_o");
    input->getItem(z_t,"length","z_t");
    input->getItem(z_m,"length","z_m");
    input->getItem(z_s,"length","z_s");
                       
    // soil section
    input->getItem(nsoilz,"soil","nsoil");
    input->getItem(soil_param,"soil","param");
    input->getItem(soil_model,"soil","model");
    input->getItem(soil_z,"soil","soil_z");
    input->getItem(soil_type,"soil","soil_type");
    input->getItem(soil_T,"soil","soil_T");
    input->getItem(soil_q,"soil","soil_q");
                       
    // initialize history arrays for first run
    soil_T_last = soil_T;
    soil_q_last = soil_q;
    
    // modify soil levels to be negative
    std::transform(soil_z.begin(), soil_z.end(), soil_z.begin(),
                   bind2nd(std::multiplies<double>(), -1.0));

    // set soil properties
    setSoilProperties();
    
    // radiation section
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
    }
                       
    // output section
    input->getItem(save_output, "output", "save");
    if (save_output) {
        
        // get fields to save from user
        input->getItem(output_fields,"output","fields");
        
        if (output_fields[0]=="all") {
            output_fields.erase(output_fields.begin());
            output_fields = {"ust","shf","lhf","ghf","obl","soilt","soilq"};
        }
        
        // create an output instance
        (j==0 && i==0) ? master=true : master=false;
        
        // add dimensions
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

        // attributes for each field
        AttScalar att_time  = {&runtime, "time",  "time",               "s",     dim_scalar_t};
        AttScalar att_ust   = {&ustar,   "ust",   "friction velocity",  "m s-1", dim_scalar_t};
        AttScalar att_shf   = {&flux_wT, "shf",   "sensible heat flux", "W m-2", dim_scalar_t};
        AttScalar att_lhf   = {&flux_wq, "lhf",   "latent heat flux",   "W m-2", dim_scalar_t};
        AttScalar att_ghf   = {&flux_gr, "ghf",   "ground heat flux",   "W m-2", dim_scalar_t};
        AttScalar att_obl   = {&L,       "obl",   "Obukhov length",     "m",     dim_scalar_t};
        AttVector att_soilz = {&soil_z,  "soilz", "Obukhov length",     "m",    dim_scalar_z};
        AttVector att_soilt = {&soil_T,  "soilt", "soil temperature",   "K",     dim_vector};
        AttVector att_soilq = {&soil_q,  "soilq", "soil moisture",      "m3 m-3",dim_vector};
        
        // map the name to attributes
        map_att_scalar.emplace("time", att_time);
        map_att_scalar.emplace("ust",  att_ust);
        map_att_scalar.emplace("shf",  att_shf);
        map_att_scalar.emplace("lhf",  att_lhf);
        map_att_scalar.emplace("ghf",  att_ghf);
        map_att_scalar.emplace("obl",  att_obl);
        map_att_vector.emplace("soilz",att_soilz);
        map_att_vector.emplace("soilt",att_soilt);
        map_att_vector.emplace("soilq",att_soilq);
        
        // we will always save time and depth
        output_scalar.push_back(map_att_scalar["time"]);
        output_vector.push_back(map_att_vector["soilz"]);
        
        // create list of fields to save
        for (int i=0; i<output_fields.size(); i++) {
            std::string key = output_fields[i];
            if (map_att_scalar.count(key)) {
                output_scalar.push_back(map_att_scalar[key]);
            } else if (map_att_vector.count(key)) {
                output_vector.push_back(map_att_vector[key]);
            }
        }
        
        // add 1D fields
        for (int i=0; i<output_scalar.size(); i++) {
            AttScalar att = output_scalar[i];
            if (master) {
                output->addField(att.name, att.units, att.long_name, att.dimensions);
            }
        }
        // add 2D fields
        for (int i=0; i<output_vector.size(); i++) {
            AttVector  att = output_vector[i];
            if (master) {
                output->addField(att.name, att.units, att.long_name, att.dimensions);
            }
        }
    }
}

// Update user-supplied fields
void UtahLSM :: updateFields(double dt,double u,double T,double q,double p,double rad=0) {
    
    tstep = dt;
    atm_U = u;
    atm_T = T;
    atm_q = q;
    atm_p = p;
    R_net = rad;
    runtime += tstep;
    
    // update total time for rad model
    if (comp_rad==1) {
        utc = std::fmod(runtime,86400);
        julian_day = int(runtime/86400);
    }
    
    // keep winds from being exactly zero
    if (atm_U==0) atm_U = 0.1;
}

void UtahLSM :: run() {
    
    surf_T_last = soil_T[0];
    surf_q_last = soil_q[0];
    
    // run radiation model is needed
    if (comp_rad==1) {
        computeRadiation();
    }
    
    // solve the surface energy balance
    solveSEB();

    // solve the surface moisture balance
    solveSMB();

    // save current temperature and moisture
    soil_T_last = soil_T;
    soil_q_last = soil_q;

    // solve diffusion equations
    solveDiffusion(1);
    solveDiffusion(2);
    
    if (first) first=false;

}

void UtahLSM :: save(Output* output) {
    
    // output size and location
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
        scalar_index = {static_cast<unsigned long>(output_counter), static_cast<unsigned long>(j)};
        vector_index = {static_cast<size_t>(output_counter), 0, static_cast<unsigned long>(j)};
        vector_size  = {1, static_cast<unsigned long>(nsoilz), 1};
        vector_index_z = {0, static_cast<unsigned long>(j)};
        vector_size_z  = {static_cast<unsigned long>(nsoilz),1};
    } else if (ny==1 && nx>1) {
        scalar_index = {static_cast<unsigned long>(output_counter), static_cast<unsigned long>(i)};
        vector_index = {static_cast<size_t>(output_counter), 0, static_cast<unsigned long>(i)};
        vector_size  = {1, static_cast<unsigned long>(nsoilz), 1};
        vector_index_z = {0, static_cast<unsigned long>(i)};
        vector_size_z  = {static_cast<unsigned long>(nsoilz),1};
    }
    else {
        scalar_index = {static_cast<unsigned long>(output_counter), static_cast<unsigned long>(j),
            static_cast<unsigned long>(i)};
        scalar_size  = {1, 1, 1};
        vector_index = {static_cast<size_t>(output_counter), 0, static_cast<unsigned long>(j),
                        static_cast<unsigned long>(i)};
        vector_size  = {1, static_cast<unsigned long>(nsoilz),1, 1};
        vector_index_z = {0, static_cast<unsigned long>(j),static_cast<unsigned long>(i)};
        vector_size_z  = {static_cast<unsigned long>(nsoilz),1 ,1};
    }
    
    // loop through 1D fields to save
    for (int i=0; i<output_scalar.size(); i++) {
        if (i==0 && master) {
            output->saveField1D(output_scalar[i].name, scalar_index, output_scalar[i].data);
        } else if (i>0) {
            output->saveField1D(output_scalar[i].name, scalar_index, output_scalar[i].data);
        }
    }
    // loop through 2D fields to save
    for (int i=0; i<output_vector.size(); i++) {
        
        // soil depth is only saved once with no time component
        if (i==0 && output_counter==0) {
            if (ny==1 && nx==1)  {
                output->saveField2D(output_vector[i].name, *output_vector[i].data);
            } else {
                output->saveField2D(output_vector[i].name, vector_index_z,
                                    vector_size_z, *output_vector[i].data);
            }
            output_vector.erase(output_vector.begin());
        } else {
            output->saveField2D(output_vector[i].name, vector_index,
                                vector_size, *output_vector[i].data);
        }
    }
    
    // increment for next time insertion
    output_counter +=1;
}

///////////////////////
// Private functions //
///////////////////////

// Set soil properties at each depth
void UtahLSM :: setSoilProperties() {
    
    struct soil::Properties soilProperties = soil::properties(soil_type,nsoilz,soil_param);
    
    b        = soilProperties.b;
    psi_sat  = soilProperties.psi_sat;
    porosity = soilProperties.porosity;
    residual = soilProperties.residual;
    K_sat    = soilProperties.K_sat;
    Ci       = soilProperties.Ci;
}

// compute net radiation
void UtahLSM :: computeRadiation() {
    R_net = radiation::net(soil_T[0],emissivity,julian_day,utc,latitude,longitude,albedo);
}

// compute surface fluxes using Monin-Obukhov w/Dyer-Hicks
void UtahLSM :: computeFluxes(double sfc_T, double sfc_q) {
    
    // local variables
    int max_iterations = 200;
    bool converged = false;
    double gnd_q, flux_wTv;
    double last_L, criteria = 0.1, ref_T = 300.;
    int depth = nsoilz;
    
    int int_depth;
    double heat_cap,dT, dz;
    struct soil::ThermalTransfer transfer;
    std::vector<double> K(depth);
    
    // compute surface mixing ratio
    gnd_q  = soil::surfaceMixingRatio(psi_sat[0],porosity[0],residual[0],b[0],
                                      sfc_T,sfc_q,atm_p,soil_model);
    
    // sensible flux, latent flux, ustar, and L
    for (int i=0; i<max_iterations; ++i) {
        
        // compute ground flux
        // first time through we estimate based on Santanello and Friedl (2003)
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
            // next time we integrate to depth of minimum flux within the soil
            transfer   = soil::thermalTransfer(psi_sat,porosity,residual,soil_q,b,Ci,depth,soil_model);
            K = transfer.k;
            
            //compute heat flux within soil and find depth of minimum
            //this is the depth to integrate time change of T
            double hfs;
            double min_sflux = 10000.;
            for (int d=0; d<nsoilz-1; ++d) {
                hfs = (d==0) ? 0.5*(K[d]+K[d+1]) * (sfc_T     - soil_T[d+1])/(soil_z[d]-soil_z[d+1]):
                0.5*(K[d]+K[d+1]) * (soil_T[d] - soil_T[d+1])/(soil_z[d]-soil_z[d+1]);
                if (std::abs(hfs)<std::abs(min_sflux)) {
                    min_sflux = hfs;
                    int_depth = d+1;
                }
            }
            flux_gr = 0;
            // Integrate time change to depth of minimum flux
            for (int d=0; d<int_depth; ++d) {
                
                heat_cap = soil::heatCapacity(porosity[d], Ci[d], soil_q[d]);
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
        
        // compute friction velocity
        ustar = atm_U*most::fm(z_m/z_o,zeta_m,zeta_o);
        
        // compute heat flux
        flux_wT = (sfc_T-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // compute latent flux
        if ( (first) && (i == 0)) {
            flux_wq = (R_net - flux_gr - flux_wT*c::rho_air*c::Cp_air)/(c::rho_air*c::Lv);
            gnd_q = atm_q + flux_wq / (ustar*most::fh(z_s/z_t,zeta_s,zeta_t));
            soil_q[0] = soil::surfaceWaterContentEstimate(psi_sat[0], porosity[0], residual[0],b[0],soil_T[0],gnd_q, atm_p,soil_model);
        } else {
            flux_wq = (gnd_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        }
        
        // compute virtual heat flux
        flux_wTv = flux_wT + ref_T*0.61*flux_wq;
        
        // compute L
        last_L = L;
        L      = -std::pow(ustar,3.)*ref_T/(c::vonk*c::grav*flux_wTv);
        
        // bounds check on L
        if (z_m/L > 5.)  L = z_m/5.;
        if (z_m/L < -5.) L = -z_m/5.;
        
        // update zeta terms
        zeta_m = z_m/L;
        zeta_s = z_s/L;
        zeta_o = z_o/L;
        zeta_t = z_t/L;
        
        // check for convergence
        converged = std::abs(last_L-L) <= criteria;
        if (converged) {
            break;
        }
    }
    
    if (!converged) {
        throw(1);
    }
}

// Solve the surface energy balance
void UtahLSM :: solveSEB() {
    
    // local variables
    int max_iter_temp = 200;
    int max_iter_flux = 200;
    double dTs, dTs_old;
    double SEB, dSEB_dT, SEB_l, SEB_h;
    double temp_l, temp_h, last_T, last_F;
    double temp_1 = soil_T[0] - 1;
    double temp_2 = soil_T[0] + 1;
    double temp_criteria = 0.001;
    double flux_criteria = 0.001;
    
    // compute SEB using current bracketed temperatures
    SEB_l = computeSEB(temp_1);
    SEB_h = computeSEB(temp_2);
    
    // dynamic bracket adjustments
    bool out_of_bracket = (SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0);
    
    while (out_of_bracket) {
        
        // expand brackets by 5 K
        temp_1 -= 1;
        temp_2 += 1;
        
        // recompute SEB at brackets
        SEB_l = computeSEB(temp_1);
        SEB_h = computeSEB(temp_2);
        
        // check for proper brackets
        out_of_bracket = (SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0);
    }
    
    if ((SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0)) {
        throw("Please adjust brackets for Ts");
    }
    
    // if SEB from low bracket Ts = 0, then that value of Ts is solution
    if (SEB_l == 0.0) soil_T[0] = temp_1;
    
    // if SEB from high bracket Ts = 0, then that value of Ts is solution
    if (SEB_h == 0.0) soil_T[0] = temp_2;
    
    // orient the solutions such that SEB(temp_l) < 0;
    if (SEB_l < 0.0) {
        temp_l = temp_1;
        temp_h = temp_2;
    } else {
        temp_l = temp_2;
        temp_h = temp_1;
    }
    
    // prepare for convergence looping
    dTs = std::abs(temp_h-temp_l);
    
    // convergence loop for flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // convergence loop for temperature
        for (int tt = 0; tt < max_iter_temp; tt++) {
            
            // compute SEB and dSEB_dTs
            SEB     = computeSEB(soil_T[0]);
            dSEB_dT = computeDSEB(soil_T[0]);
            
            // update brackets
            if (SEB<0.) temp_l = soil_T[0];
            if (SEB>0.) temp_h = soil_T[0];
            
            // bracket and bisect temperature if Newton out of range
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
            
            // check for convergence
            if (std::abs( (soil_T[0]-last_T)/last_T) <= temp_criteria) break;
            
            // if convergence fails, recompute flux
            computeFluxes(soil_T[0],soil_q[0]);
        }
        
        // save current flux for convergence criteria
        last_F = flux_wT;
        
        // recompute heat flux using new temperature
        computeFluxes(soil_T[0],soil_q[0]);
        
        // check for convergence
        if (std::abs(flux_wT-last_F) <= flux_criteria) {
            break;
        }
        
        // if flux fails to converge, split temperature
        soil_T[0] = 0.5*(soil_T[0] + last_T);
        
        // recompute stability corrections and fluxes
        computeFluxes(soil_T[0],soil_q[0]);
    }
}

// computes SEB and dSEB_dTs
double UtahLSM :: computeSEB(double sfc_T) {
    
    // local variables
    int depth = nsoilz;
    double Qg, Qh, Ql, SEB;
    std::vector<double> K_soil(depth);
    
    // compute fluxes using passed in values
    computeFluxes(sfc_T,soil_q[0]);
    
    // write sensible and latent heat fluxes in [W/m^2]
    Qh = c::rho_air*c::Cp_air*flux_wT;
    Ql = c::rho_air*c::Lv*flux_wq;
    Qg = flux_gr;

    // compute surface energy balance
    SEB = R_net - Qg - Qh - Ql;
    return SEB;
}

// computes derivative of SEB wrt surface temperature
double UtahLSM :: computeDSEB(double sfc_T) {
    
    // local variables
    double heat_cap, dSEB_dT;
    
    //compute derivative of SEB wrt temperature
    heat_cap = soil::heatCapacity(porosity[0], Ci[0], soil_q[0]);
    dSEB_dT = - 4.*emissivity*c::sb*std::pow(sfc_T,3.)
    - c::rho_air*c::Cp_air*ustar*most::fh(z_s/z_t,zeta_s,zeta_t)
    - heat_cap*(soil_z[0]-soil_z[1])/(2*tstep);
    return dSEB_dT;
}

// Solve the surface energy balance
void UtahLSM :: solveSMB() {
    
    // local variables
    bool converged;
    int max_iter_flux = 200;
    double E;
    double flux_sm_last, flux_sm, flux_sm2;
    double K_n_avg, D_n_avg;
    double delta = 0.8, flux_criteria = .001;
    std::vector<double> psi;
    std::vector<double> D_n;
    std::vector<double> K_n;
    struct soil::MoistureTransfer transfer;
    
    // moisture potential at first level below ground
    psi = soil::waterPotential(psi_sat, porosity, residual, soil_q, b, nsoilz, soil_model);
    
    // compute initial soil moisture flux
    transfer = soil::moistureTransfer(psi_sat,K_sat,porosity,residual,soil_q,b,2,soil_model);
    K_n      = transfer.k;
    D_n      = transfer.d;
    K_n_avg  = std::accumulate(K_n.begin(), K_n.end(), 0.0)/K_n.size();
    D_n_avg  = std::accumulate(D_n.begin(), D_n.end(), 0.0)/D_n.size();
    flux_sm  = c::rho_wat*K_n_avg*((psi[0] - psi[1])/(soil_z[0]-soil_z[1]) + 1);
    flux_sm2 = c::rho_wat*D_n_avg*(soil_q[0]-soil_q[1])/(soil_z[0]-soil_z[1])
    + c::rho_wat*K_n_avg;
    
    // compute evaporation
    E = c::rho_air*flux_wq;
    
    // convergence loop for moisture flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // save soil moisture flux for convergence test
        flux_sm_last = flux_sm;
        
        // compute new weighted soil moisture flux
        flux_sm = delta*flux_sm_last - (1.-delta)*E;
        
        // re-compute moisture potential
        psi[0]    = psi[1] + (soil_z[0]-soil_z[1])*((flux_sm/(c::rho_wat*K_n_avg))-1);
        if (psi[0]>psi_sat[0]) {
            psi[0] = psi_sat[0];
        }
        
        // update soil moisture
        soil_q[0] = soil::surfaceWaterContent(psi[0], psi_sat[0], porosity[0],
                                              residual[0], b[0], soil_model);

        double gnd_q  = soil::surfaceMixingRatio(psi_sat[0],porosity[0],residual[0],b[0],
                                                 soil_T[0],soil_q[0],atm_p,soil_model);
        E = c::rho_air*(gnd_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // update soil moisture transfer
        transfer = soil::moistureTransfer(psi_sat,K_sat,porosity,residual,soil_q,b,2,soil_model);
        K_n      = transfer.k;
        K_n_avg  = std::accumulate(K_n.begin(), K_n.end(), 0.0)/K_n.size();
        
        // check for convergence
        converged = std::abs((E + flux_sm)/E) <=flux_criteria;
    }
}

// integrate soil heat diffusion equation
void UtahLSM :: solveDiffusion(int type) {
    
    // local variables
    double surf_scalar_last;
    double dKdz, C_tp, C_tm, C_p, C_m;
    double dz_p, dz_m, a_m, a_p, a_o;
    struct soil::MoistureTransfer transfer_m;
    struct soil::ThermalTransfer transfer_t;
    std::vector<double> K(nsoilz,0.0);
    std::vector<double> D(nsoilz,0.0);
    std::vector<double> K_mid(nsoilz-1,0.0);
    std::vector<double> z_mid(nsoilz-1,0.0);
    std::vector<double> r(nsoilz-1,0.0);
    std::vector<double> e(nsoilz-1,0.0);
    std::vector<double> f(nsoilz-1,0.0);
    std::vector<double> g(nsoilz-1,0.0);
    std::vector<double> scalar(nsoilz);
    
    // set up and solve a tridiagonal matrix
    // AT(n+1) = r(n), where n denotes the time level
    // e, f, g the components of A matrix
    // T(n+1)  the soil temperature vector at t=n+1
    // r(n)    the soil temperature vector at t=n multiplied by coefficients
    
    // interpolate soil_z and D_n to mid-points
    if (type==1) {
        transfer_t = soil::thermalTransfer(psi_sat,porosity,porosity,soil_q,b,Ci,nsoilz,soil_model);
        D = transfer_t.d;
        for (int i=0; i<nsoilz-1; i++) {
            K_mid[i] = 0.5*(D[i]+D[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar = soil_T;
        surf_scalar_last = soil_T_last[0];
    } else {
        transfer_m = soil::moistureTransfer(psi_sat,K_sat,porosity,residual,soil_q,b,nsoilz,soil_model);
        D      = transfer_m.d;
        K      = transfer_m.k;
        for (int i=0; i<nsoilz-1; i++) {
            K_mid[i] = 0.5*(D[i]+D[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar = soil_q;
        surf_scalar_last = soil_q_last[0];
    }
    
    // matrix coefficients for first level below surface
    C_p  = tstep*K_mid[0]/(2*(z_mid[0]-z_mid[1])*(soil_z[0]-soil_z[1]));
    C_m  = tstep*K_mid[1]/(2*(z_mid[0]-z_mid[1])*(soil_z[1]-soil_z[2]));
    C_tp = (1.0 + C_p + C_m);
    C_tm = (1.0 - C_p - C_m);
    
    f[0] =  C_tp;
    g[0] = -C_m;
    r[0] =  C_p*surf_scalar_last + C_tm*scalar[1] + C_m*scalar[2] + C_p*scalar[0];
    
    // for moisture, we need additional dKn/dz term
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
    
    // matrix coefficients for the interior levels
    for (int i=1; i<nsoilz-2; i++) {
        C_p  = tstep*K_mid[i]  / (2*(z_mid[i]-z_mid[i+1])*(soil_z[i]  -soil_z[i+1]));
        C_m  = tstep*K_mid[i+1]/ (2*(z_mid[i]-z_mid[i+1])*(soil_z[i+1]-soil_z[i+2]));
        C_tp = (1 + C_p + C_m);
        C_tm = (1 - C_p - C_m);
        
        e[i] = -C_p;
        f[i] =  C_tp;
        g[i] = -C_m;
        r[i] =  C_p*scalar[i] + C_tm*scalar[i+1] + C_m*scalar[i+2];
        
        // for moisture, we need additional dKn/dz term
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
    
    // matrix coefficients for bottom level
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
    
    // for moisture, we need additional dKn/dz term
    // use simple 2-pt backward Euler approximation
    if (type==2) {
        dKdz = tstep*(K[nsoilz-2]-K[nsoilz-1])/(soil_z[nsoilz-2]-soil_z[nsoilz-1]);
        r[nsoilz-2] = r[nsoilz-2] + dKdz;
    }
    
    // solve the tridiagonal system
    try {
        if (type==1) matrix::tridiagonal(e,f,g,r,soil_T);
        if (type==2) matrix::tridiagonal(e,f,g,r,soil_q);
    } catch(std::string &e) {
        std::cout<<e<<std::endl;
        std::exit(0);
    }
}
