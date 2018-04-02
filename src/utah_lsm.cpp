//
//  utah_lsm.cpp
//  
//  This class handles the Utah LSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>
#include <functional>
#include "utah_lsm.hpp"
#include "constants.hpp"
#include "soil.hpp"
#include "radiation.hpp"
#include "most.hpp"
#include "matrix.hpp"

namespace {
    namespace c = Constants;
}

UtahLSM::UtahLSM(bool first, double dt, double z_o,double z_t,double z_m,double z_s, 
                 double atm_p,double atm_ws,double atm_T,double atm_q, 
                 int nsoilz,std::vector<double>& soil_z,std::vector<int>& soil_type,
                 std::vector<double>& soil_T,std::vector<double>& soil_q,
                 int julian_day, double utc, double latitude, double longitude,
                 double albedo, double emissivity, double R_net, int comp_rad,
                 double& zeta_m,double& zeta_s,double& zeta_o,double& zeta_t,
                 double& ustar, double& flux_wT, double& flux_wq) : 
                 first(first),dt(dt),z_o(z_o),z_t(z_t),z_m(z_m),z_s(z_s),
                 atm_p(atm_p),atm_ws(atm_ws),atm_T(atm_T),atm_q(atm_q),
                 nsoilz(nsoilz),soil_z(soil_z),soil_type(soil_type),soil_T(soil_T),soil_q(soil_q),
                 julian_day(julian_day),utc(utc),latitude(latitude),longitude(longitude), 
                 albedo(albedo),emissivity(emissivity),R_net(R_net),comp_rad(comp_rad),
                 zeta_m(zeta_m),zeta_s(zeta_s),zeta_o(zeta_o),zeta_t(zeta_t),
                 ustar(ustar),flux_wT(flux_wT),flux_wq(flux_wq) {
        
    // local variables
    surf_T_last = soil_T[0];
    surf_q_last = soil_q[0];
    
    // run the model
    defineSoil();
    if (first) computeFluxes(3); 
    if (comp_rad==true) computeRadiation();
    solveSEB();
    solveMoisture();
    solveDiffusion(1);
    solveDiffusion(2);
}

// Construct soil properties
void UtahLSM :: defineSoil() {
	
	soil::soilProperties soilProperties = soil::soilTypeProperties(soil_type,nsoilz);
	b        = soilProperties.b;
	psi_sat  = soilProperties.psi_sat;
	porosity = soilProperties.porosity;
	K_sat    = soilProperties.K_sat;
	Ci       = soilProperties.Ci;
}

// compute surface fluxes using Monin-Obukhov w/Dyer-Hicks
// flux = 1: update sensible heat flux
// flux = 2: update latent heat flux
// flux = 3: update both fluxes
void UtahLSM :: computeFluxes(int flux) {
    
    // local variables
    int max_iterations = 200;
    bool converged = false;
    double L=0, q_gnd, q_zot, T_gnd, T_zot,flux_wTv, atm_Tv;
    double last_L, criteria = 0.1, ref_T = 300.;
    
    // compute surface mixing ratio
    q_gnd = soil::surfaceMixingRatio(psi_sat[0],porosity[0],b[0],
                                     soil_T[0],soil_q[0],atm_p);
    
    //std::cout<<soil_q[0]<<std::endl;
    
    atm_Tv = atm_T * (1 + 0.61*atm_q);
    
    //zeta_s = 0, zeta_t = 0, zeta_o = 0, zeta_m = 0;
    
    // iterate for convergence
    for (int i=0; i<max_iterations; ++i) {
        
        // approximate ground temperature/moisture to roughness height
        T_zot = soil_T[0] + (atm_T-soil_T[0])*std::log(z_o/z_t)*most::fh(z_s/z_t,zeta_s,zeta_t)/c::vonk;
        q_zot = q_gnd + (atm_q-q_gnd)*std::log(z_o/z_t)*most::fh(z_s/z_t,zeta_s,zeta_t)/c::vonk;
         
        //T_zot = soil_T[0];
        //q_zot = q_gnd;
              
        // compute friction velocity
        ustar = atm_ws*most::fm(z_m/z_o,zeta_m,zeta_o);
         
        // compute heat flux
        if ( flux==1 || flux==3 ){
            flux_wT = (T_zot-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        }
        
        // compute latent heat flux
        if ( flux==2 || flux==3 ){
            flux_wq = (q_zot-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        }
        
        // compute virtual heat flux
        flux_wTv = flux_wT + ref_T*0.61*flux_wq;
                
        // compute L
        last_L = L;
        L = -std::pow(ustar,3.)*ref_T/(c::vonk*c::grav*flux_wTv);
        
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
        if (converged) break;
    }
    
    if (!converged) {
	    std::cout<<"Something is wrong: L didn't converge"<<std::endl;
	    throw(1);
	}
}

// compute net radiation using simple model
void UtahLSM :: computeRadiation() {

    R_net = radiation::net(soil_T[0],emissivity,julian_day,utc,
                           latitude,longitude,albedo);
}

// Solve the surface energy balance
void UtahLSM :: solveSEB() {
        
    // local variables
    int max_iter_temp = 150;
    int max_iter_flux = 150; 
    double dTs, dTs_old;
    double SEB, dSEB_dT, SEB_l, SEB_h;
    double temp_l, temp_h, last_T, last_F;
    double temp_1 = soil_T[0] - 150;
    double temp_2 = soil_T[0] + 150;
    double temp_criteria = 0.001;
    double flux_criteria = 0.0001;
        
    // compute SEB using current+bracketed temperatures
    SEB_l = computeSEB(temp_1);
    SEB_h = computeSEB(temp_2);
    
    // check for improper bracketing
    if ((SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0)) {
        throw("Please adjust brackets for Ts");
    }
    
    // if SEB from low bracket Ts = 0, then that value of Ts is solution
    if (SEB_l == 0.0) soil_T[0] = temp_1;

    // if SEB from high bracket Ts = 0, then that value of Ts is solution
    if (SEB_h == 0.0) soil_T[0] = temp_2;
    
    // orient the solutions, such that SEB(temp_l) < 0;
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
            if (std::abs(dTs) <= temp_criteria) break;
            
            // if convergence fails, recompute flux
            computeFluxes(3);
        }
        
        // save current flux for convergence criteria
        last_F = flux_wT;
        
        // recompute heat flux using new temperature
        computeFluxes(3);
        
        // check for convergence
        if (std::abs(flux_wT-last_F) <= flux_criteria) break;
                
        // if flux fails to converge, split temperature
        soil_T[0] = 0.5*(soil_T[0] + last_T);
        
        // recompute stability corrections and fluxes
        computeFluxes(3);
    }
}

// computes SEB and dSEB_dTs
double UtahLSM :: computeSEB(double sfc_T) {
    
    // local variables
    int depth = 2;
    double K_soil_avg, Qg, Qh, Ql, SEB, C;
    double dQh_dT, dSEB_dT, dTs, Rl_up;
    std::vector<double> K_soil(depth);
    
    // compute soil heat flux
    K_soil     = soil::soilThermalTransfer(psi_sat,porosity,soil_q,b,Ci,depth,0);
    K_soil_avg = std::accumulate(K_soil.begin(), K_soil.end(), 0.0)/K_soil.size();
    Qg         = K_soil_avg*(sfc_T - soil_T[1])/(soil_z[0]-soil_z[1]);
    
    // write sensible and latent heat fluxes in [W/m^2]
    Qh = c::rho_air*c::Cp_air*flux_wT;
    Ql = c::rho_air*c::Lv*flux_wq;
    
    // compute surface energy balance
    SEB = R_net - Qg - Qh - Ql;
    
    return SEB;
}

// computes derivative of SEB wrt surface temperature
double UtahLSM :: computeDSEB(double sfc_T) {
    
    // local variables
    int depth = 2;
    double K_soil_avg, dSEB_dT;
    std::vector<double> K_soil;
    
    // compute soil heat transfer
    K_soil     = soil::soilThermalTransfer(psi_sat,porosity,soil_q,b,Ci,depth,0);
    K_soil_avg = std::accumulate(K_soil.begin(), K_soil.end(), 0.0)/K_soil.size();
    
    //compute derivative of SEB wrt temperature
    dSEB_dT = - 4.*emissivity*c::sb*std::pow(sfc_T,3.)
              - K_soil_avg/(soil_z[0]-soil_z[1]) 
              - c::rho_air*c::Cp_air*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);;
    
    return dSEB_dT;
}

// solve for moisture flux
void UtahLSM :: solveMoisture() {
    
    // local variables
    int max_iter_flux = 150;
    double flux_sm_last, flux_sm, flux_sm2;
    double psi_n0, psi_n1, sfc_q, K_n_avg, D_n_avg;
    double delta = 0.85, flux_criteria = .000001;
    std::vector<double> D_n;
    std::vector<double> K_n;
    soil::soilTransfer transfer;
    
    // moisture potential at first level below ground
    psi_n0 = psi_sat[0]*std::pow(porosity[0]/soil_q[0],b[0]);
    psi_n1 = psi_sat[1]*std::pow(porosity[1]/soil_q[1],b[1]);
    
    // compute initial soil moisture flux
    transfer = soil::soilMoistureTransfer(psi_sat,K_sat,porosity,soil_q,b,2);
    K_n      = transfer.transfer_h;
    D_n      = transfer.transfer_d;
    K_n_avg  = std::accumulate(K_n.begin(), K_n.end(), 0.0)/K_n.size();
    D_n_avg  = std::accumulate(D_n.begin(), D_n.end(), 0.0)/D_n.size();
    flux_sm  = c::rho_wat*K_n_avg*((psi_n0 - psi_n1)/(soil_z[0]-soil_z[1]) + 1);
    flux_sm2 = c::rho_wat*D_n_avg*(soil_q[0]-soil_q[1])/(soil_z[0]-soil_z[1])
                     + c::rho_wat*K_n_avg;
    
    std::cout<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"Soil Flux K: "<<flux_sm<<std::endl;
    std::cout<<"Soil Flux D: "<<flux_sm2<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"Psi_s  Poros Soilq1 Soilq2 Exp"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<psi_sat[0]<<" "<<porosity[0]<<" "<<soil_q[0]<<" "<<soil_q[1]<<" "<<b[0]<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    //std::cout<<psi_sat[0]<<" "<<psi_n0<<" "<<psi_n1<<" "<<soil_q[0]<<" "<<soil_q[1]<<std::endl;
    //std::cout<<K_n_avg<<" "<<flux_sm<<" "<<(psi_n0-psi_n1)/(soil_z[0]-soil_z[1])<<std::endl;    
    
    // convergence loop for moisture flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // recompute moisture flux
        computeFluxes(2);
        
        // save soil moisture flux for convergence test
        flux_sm_last = flux_sm;

        // compute new weighted flux
        flux_sm = delta*flux_sm_last - (1.-delta)*c::rho_air*flux_wq;
        
        std::cout<<"Old soil flux: "<<flux_sm_last<<std::endl;
        std::cout<<"New soil flux: "<<flux_sm<<std::endl;
        std::cout<<"---------------------------------"<<std::endl;
        // update surface moisture
        std::cout<<"Old psi_n0: "<<psi_n0<<std::endl;
        psi_n0    = psi_n1 + (soil_z[0]-soil_z[1])*((flux_sm/(c::rho_wat*K_n_avg))-1);
        std::cout<<"New psi_n0: "<<psi_n0<<std::endl;
        std::cout<<"---------------------------------"<<std::endl;
        std::cout<<"Psi_1    sz0 sz1  Kbar"<<std::endl;
		std::cout<<"---------------------------------"<<std::endl;
        std::cout<<psi_n1<<" "<<soil_z[0]<<" "<<soil_z[1]<<" "<<K_n_avg<<std::endl;
        std::cout<<"---------------------------------"<<std::endl;
        std::cout<<"Old soil_q[0]: "<<soil_q[0]<<std::endl;
        soil_q[0] = porosity[0]*std::pow(psi_sat[0]/psi_n0,(1./b[0]));
        std::cout<<"New soil_q[0]: "<<soil_q[0]<<std::endl;
        std::cout<<"---------------------------------"<<std::endl;
        std::cout<<"Psi_s  Psi_0   Poros Soilq1 Soilq2 Exp"<<std::endl;
		std::cout<<"---------------------------------"<<std::endl;
		std::cout<<psi_sat[0]<<" "<<psi_n0<<" "<<porosity[0]<<" "<<soil_q[0]<<" "<<soil_q[1]<<" "<<b[0]<<std::endl;        
        std::cout<<"---------------------------------"<<std::endl;
        
        // update soil moisture transfer
        transfer = soil::soilMoistureTransfer(psi_sat,K_sat,porosity,soil_q,b,2);
        K_n      = transfer.transfer_h;
        K_n_avg  = std::accumulate(K_n.begin(), K_n.end(), 0.0)/K_n.size(); 
        
        // check for convergence
        if (std::abs((-c::rho_air*flux_wq - flux_sm)
           / (-c::rho_air*flux_wq)) <=flux_criteria) {
             break;  
        }           
    }
    std::cout<<"--- Converged with heat flux = "<<flux_wT<<", mois flux = "<<flux_wq<<", and temp = "<<soil_T[0]<<std::endl;
}

// integrate soil heat diffusion equation
void UtahLSM :: solveDiffusion(int type) {
    
    // local variables
    double surf_scalar_last;
    double dKdz, C_tp, C_tm, C_p, C_m;
    double dz_p, dz_m, a_m, a_p, a_o;
    soil::soilTransfer transfer;
    std::vector<double> K_n(nsoilz,0.0);
    std::vector<double> D_n(nsoilz,0.0);
    std::vector<double> K_mid(nsoilz-1,0.0);
    std::vector<double> z_mid(nsoilz-1,0.0);
    std::vector<double> r(nsoilz-1,0.0);
    std::vector<double> e(nsoilz-1,0.0);
    std::vector<double> f(nsoilz-1,0.0);
    std::vector<double> g(nsoilz-1,0.0);
    std::vector<double> scalar(nsoilz);
    
    // here we set up and solve a tridiagonal matrix
    // AT(n+1) = r(n), where n denotes the time level
    // e, f, g the components of A matrix
    // T(n+1)  the soil temperature vector at t=n+1
    // r(n)    the soil temperature vector at t=n multiplied by coefficients         
    
    // interpolate soil_z and D_n to mid-points 
    if (type==1) {
        D_n = soil::soilThermalTransfer(psi_sat,porosity,soil_q,b,Ci,nsoilz,1);
        for (int i=0; i<nsoilz-1; i++) {
            K_mid[i] = 0.5*(D_n[i]+D_n[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar = soil_T;
        surf_scalar_last = surf_T_last;
    } else {
        transfer = soil::soilMoistureTransfer(psi_sat,K_sat,porosity,soil_q,b,nsoilz);
        D_n      = transfer.transfer_d;
        K_n      = transfer.transfer_h;
        for (int i=0; i<nsoilz-1; i++) {
            K_mid[i] = 0.5*(D_n[i]+D_n[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar = soil_q;
        surf_scalar_last = surf_q_last;
    }
    
    // matrix coefficients for first level below surface
    C_p  = dt*K_mid[0]/(2*(z_mid[0]-z_mid[1])*(soil_z[0]-soil_z[1]));
    C_m  = dt*K_mid[1]/(2*(z_mid[0]-z_mid[1])*(soil_z[1]-soil_z[2]));
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
        
        dKdz = dt*(K_n[0]*a_p + K_n[1]*a_o + K_n[2]*a_m);
        r[0] = r[0] + dKdz;          
    }
    
    // matrix coefficients for the interior levels        
    for (int i=1; i<nsoilz-2; i++) {
        C_p  = dt*K_mid[i]  / (2*(z_mid[i]-z_mid[i+1])*(soil_z[i]  -soil_z[i+1]));
        C_m  = dt*K_mid[i+1]/ (2*(z_mid[i]-z_mid[i+1])*(soil_z[i+1]-soil_z[i+2]));
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
            dKdz = dt*(K_n[i]*a_p + K_n[i+1]*a_o + K_n[i+2]*a_m);
            r[i] = r[i] + dKdz;
        }
    }
    
    // matrix coefficients for bottom level
    // first, construct ghost values
    int j       = nsoilz-2;
    double z_g  = 2*soil_z[j+1] - soil_z[j];
    double z_mg = (soil_z[j+1] + z_g) / 2.;
    
    C_p  = dt*K_mid[j]/(2*(z_mid[j]-z_mg)*(soil_z[j]-soil_z[j+1]));
    C_m  = dt*K_mid[j]/(2*(z_mid[j]-z_mg)*(soil_z[j+1]-z_g));
    C_tp = (1 + C_p + C_m);
    C_tm = (1 - C_p - C_m);
    
    e[j] = C_m  - C_p;
    f[j] = C_tp - 2* C_m;
    r[j] = (C_p - C_m)*scalar[j] + (C_tm+2*C_m)*scalar[j+1];
    
    // for moisture, we need additional dKn/dz term
    // use simple 2-pt backward Euler approximation 
    if (type==2) {
        dKdz = dt*(K_n[nsoilz-2]-K_n[nsoilz-1])/(soil_z[nsoilz-2]-soil_z[nsoilz-1]);
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