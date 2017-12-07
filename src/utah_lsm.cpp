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
#include "most.hpp"
#include "matrix.hpp"

namespace {
    namespace c = Constants;
}

UtahLSM::UtahLSM(double dt, double z_o,double z_t,double z_m,double z_s, 
                double atm_p,double atm_ws,double atm_T,double atm_q, 
                int nsoilz,std::vector<double>& soil_z,std::vector<double>& soil_T,
                std::vector<double>& soil_q,std::vector<double>& porosity,
                std::vector<double>& psi_nsat,std::vector<double>& K_nsat,
                std::vector<double>& b,std::vector<double>& Ci,
                int julian_day, double utc, double latitude, double longitude,
                double albedo, double emissivity, double R_net,
                double& zeta_m,double& zeta_s,double& zeta_o,double& zeta_t,
                double& ustar, double& flux_wT, double& flux_wq) : 
                dt(dt),z_o(z_o),z_t(z_t),z_m(z_m),z_s(z_s),
                atm_p(atm_p),atm_ws(atm_ws),atm_T(atm_T),atm_q(atm_q),
                nsoilz(nsoilz),julian_day(nsoilz),utc(utc),latitude(latitude), 
                longitude(longitude), albedo(albedo),emissivity(emissivity),R_net(R_net),
                soil_z(soil_z),soil_T(soil_T),soil_q(soil_q),
                porosity(porosity),psi_nsat(psi_nsat), K_nsat(K_nsat),b(b),Ci(Ci),
                zeta_m(zeta_m),zeta_s(zeta_s),zeta_o(zeta_o),zeta_t(zeta_t),
                ustar(ustar),flux_wT(flux_wT),flux_wq(flux_wq) {
        
    // local variables
    surf_T_last = soil_T[0];
    surf_q_last = soil_q[0];
    
    // run the model    
    computeFluxes();
    solveSEB();
    solveMoisture();
    solveDiffusion(1);
    solveDiffusion(2);
}

// compute surface fluxes using Monin-Obukhov w/Dyer-Hicks
void UtahLSM :: computeFluxes() {
    
    // local variables
    double L, sfc_q, flux_wTv;
    
    // compute surface mixing ratio
    sfc_q = soil::surfaceMixingRatio(psi_nsat[0],porosity[0],b[0],
                                     soil_T[0],soil_q[0],atm_p);
    
    // iterate for convergence
    for (int i=0; i<4; ++i) {
        
        // compute friction velocity
        ustar = atm_ws*most::fm(z_m/z_o,zeta_m,zeta_o);
        
        // compute heat flux
        flux_wT = (soil_T[0]-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);

        // compute latent heat flux
        flux_wq = (sfc_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // compute virtual heat flux
        flux_wTv = flux_wT + atm_T*0.61*flux_wq;
        
        // compute L
        L = std::pow(-ustar,3.)*atm_T/(c::vonk*c::grav*flux_wTv);
        
        // bounds check on L
        if (z_m/L > 5.)  L = z_m/5.;
        if (z_m/L < -5.) L = -z_m/5.;
        
        // update zeta terms
        zeta_m = z_m/L;
        zeta_s = z_s/L;
        zeta_o = z_o/L;
        zeta_t = z_t/L; 
    }
}

// Solve the surface energy balance
void UtahLSM :: solveSEB() {
    
    // local variables
    int max_iter_temp = 100;
    int max_iter_flux = 100; 
    double dTs, dTs_old;
    double SEB, dSEB_dT, temp_l, temp_h;
    double SEB_l, SEB_h, last_T, last_F;
    double temp_1 = soil_T[0] - 30;
    double temp_2 = soil_T[0] + 30;
    double temp_criteria = 0.001;
    double flux_criteria = 0.001;
        
    // compute SEB with bracketed temperatures
    SEB_l = computeSEB(temp_1);
    SEB_h = computeSEB(temp_2);
    double SEB_o = computeSEB(soil_T[0]);
    
    //std::cout<<SEB_l<<" "<<SEB_o<<" "<<SEB_h<<std::endl;
    
    // check for improper bracketing
    if ((SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0)) {
        throw("Please adjust brackets for Ts");
    }
    
    // if SEB from low bracket Ts = 0, then that value of Ts is solution
    if (SEB_l == 0.0) {
        soil_T[0] = temp_1;
    }
    // if SEB from high bracket Ts = 0, then that value of Ts is solution
    if (SEB_h == 0.0){
        soil_T[0] = temp_2;
    }
    // orient the solutions, such that SEB(temp_l) < 0;
    if (SEB_l < 0.0) {
        temp_l = temp_1;
        temp_h = temp_2;
    } else {
        temp_h = temp_1;
        temp_l = temp_2; 
    }
    
    // prepare for convergence looping
    dTs     = std::abs(temp_2-temp_1);
    dTs_old = dTs;
    SEB     = computeSEB(soil_T[0]);
    dSEB_dT = computeDSEB(soil_T[0]);
     
    // convergence loop for flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // convergence loop for temperature
        for (int tt = 0; tt < max_iter_temp; tt++) {
            
            // bracket and bisect temperature if Newton out of range
            if ((((soil_T[0]-temp_h)*dSEB_dT-SEB)*((soil_T[0]-temp_l)*dSEB_dT-SEB)>0.0) 
               || (std::abs(2*SEB) > std::abs(dTs_old*dSEB_dT))) {
                dTs_old   = dTs;
                dTs       = 0.5*(temp_h-temp_l);
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
            if (std::abs(dTs) <= temp_criteria) {
                std::cout<<"Found Temp!"<<std::endl;
                break;
            }
            
            // if convergence fails, recompute flux
            flux_wT = (soil_T[0]-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
            
            // recompute SEB and dSEB_dT since f(Ts)
            SEB     = computeSEB(soil_T[0]);
            dSEB_dT = computeDSEB(soil_T[0]);
            
            // update brackets
            if (SEB<0.) {
                temp_l = soil_T[0];
            } else {
                temp_h = soil_T[0];
            }
        }
        
        // save current flux for convergence criteria
        last_F = flux_wT;
        
        // recompute heat flux using new temperature
        flux_wT = (soil_T[0]-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
    
        // check for convergence
        if (std::abs(flux_wT-last_F) <= flux_criteria) break;
                
        // if flux fails to converge, split temperature
        soil_T[0] = 0.5*(soil_T[0] + last_T);
        last_T = soil_T[0];
        
        // recompute SEB, dSEB_dT with new temperature
        SEB     = computeSEB(soil_T[0]);
        dSEB_dT = computeDSEB(soil_T[0]);
        
        // adjust brackets
        if (SEB<0.) {
            temp_l = soil_T[0];
        } else {
            temp_h = soil_T[0];
        }
        //computeFluxes();
    }
    std::cout<<"--- Converged with heat flux = "<<flux_wT<<" and temp = "<<soil_T[0]<<std::endl;
}

// computes SEB and dSEB_dTs
double UtahLSM :: computeSEB(double sfc_T) {
    
    // local variables
    double K_soil_avg, Qg, Qh, Ql, SEB, C;
    double dQh_dT, dSEB_dT, dTs, Rl_up;
    std::vector<double> K_soil;
    
    // compute soil heat flux
    K_soil     = soil::soilThermalTransfer(psi_nsat,porosity,soil_q,b,Ci,2,0);
    K_soil_avg = std::accumulate(K_soil.begin(), K_soil.end(), 0.0)/K_soil.size();
    Qg         = K_soil_avg*(soil_T[1]-sfc_T)/(soil_z[1]-soil_z[0]);
    
    // write sensible and latent heat fluxes in [W/m^2]
    Qh = Constants::rho_air*c::Cp_air*flux_wT;
    Ql = Constants::rho_air*c::Lv*flux_wq;
    
    // Rnet is from measurements, separate lw up
    Rl_up = emissivity*c::sb*std::pow(sfc_T,4.);
    
    // add together terms that are not f(Ts)
    C = Ql -(R_net + Rl_up);
    
    // compute surface energy balance
    SEB = Qg + Qh + Rl_up + C;
    
    return SEB;
}

// computes derivative of SEB wrt surface temperature
double UtahLSM :: computeDSEB(double sfc_T) {
    
    // local variables
    double K_soil_avg, dQh_dT, dSEB_dT;
    std::vector<double> K_soil;
    
    // compute soil heat transfer
    K_soil     = soil::soilThermalTransfer(psi_nsat,porosity,soil_q,b,Ci,2,0);
    K_soil_avg = std::accumulate(K_soil.begin(), K_soil.end(), 0.0)/K_soil.size();
    
    // compute derivative of heat flux wrt temperature
    dQh_dT = c::rho_air*c::Cp_air*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
    dQh_dT = c::rho_air*c::Cp_air*dQh_dT;
    
    //compute derivative of SEB wrt temperature
    dSEB_dT = 4.*emissivity*c::sb*std::pow(sfc_T,3.)
              -K_soil_avg/(soil_z[1]-soil_z[0]) + dQh_dT;
    
    return dSEB_dT;
}

// solve for moisture flux
void UtahLSM :: solveMoisture() {
    
    // local variables
    bool fluxConverged = false;
    int max_iter_flux = 100;
    double flux_sm_last, flux_sm;
    double psi_n0, psi_n1, sfc_q;
    double D_n_avg, K_n_avg;
    double delta = .7, flux_criteria = 0.001;
    std::vector<double> D_n;
    std::vector<double> K_n;
    
    // moisture potential at first level below ground
    psi_n0 = psi_nsat[0]*std::pow(porosity[0]/soil_q[0],b[0]);
    psi_n1 = psi_nsat[1]*std::pow(porosity[1]/soil_q[1],b[1]);
    
    // compute initial soil moisture flux
    std::tie(D_n, K_n) = soil::soilMoistureTransfer(psi_nsat,K_nsat,porosity,soil_q,b,2);
    D_n_avg = std::accumulate(D_n.begin(), D_n.end(), 0.0)/D_n.size();
    K_n_avg = std::accumulate(K_n.begin(), K_n.end(), 0.0)/K_n.size();
    
    //flux_sm = -c::rho_wat*D_n*(soil_q[0]-soil_q[1])/(soil_z[1]-soil_z[0]) 
                     //+ c::rho_wat*K_n;
    flux_sm = c::rho_wat*K_n_avg*((psi_n0 - psi_n1)/(soil_z[0]-soil_z[1]) + 1);
        
    // convergence loop for moisture flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // compute surface mixing ratio
        sfc_q = soil::surfaceMixingRatio(psi_nsat[0],porosity[0],b[0],
                                         soil_T[0],soil_q[0],atm_p);
                       
        // compute atm moisture flux
        flux_wq = (sfc_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // check for convergence
        if (std::abs((c::rho_air*flux_wq + flux_sm)
           / (c::rho_air*flux_wq)) <=flux_criteria) {
             break;  
        }
        
        // save soil moisture flux for convergence test
        flux_sm_last = flux_sm;
        
        // update soil moisture flux
        flux_sm = delta*flux_sm_last - (1.-delta)*c::rho_air*flux_wq;

        // update soil moisture transfer
        std::tie(D_n, K_n) = soil::soilMoistureTransfer(psi_nsat,K_nsat,porosity,soil_q,b,2);
        D_n_avg = std::accumulate(D_n.begin(), D_n.end(), 0.0)/D_n.size();
        K_n_avg = std::accumulate(K_n.begin(), K_n.end(), 0.0)/K_n.size();
          
        // update surface moisture
        psi_n0    = psi_n1 + (soil_z[0]-soil_z[1])*((flux_sm/(c::rho_wat*K_n_avg))-1);
        soil_q[0] = porosity[0]*std::pow(psi_n0/psi_nsat[0],-(1./b[0]));                   
    }
    std::cout<<"--- Converged with mois flux = "<<flux_wq<<std::endl;
}

// integrate soil heat diffusion equation
void UtahLSM :: solveDiffusion(int type) {
    
    // local variables
    double surf_scalar_last;
    double dKdz, C_tp, C_tm, C_p, C_m;
    double dz_p, dz_m, a_m, a_p, a_o;
    std::vector<double> D_n(nsoilz,0.0);
    std::vector<double> K_n(nsoilz,0.0);
    std::vector<double> K_mid(nsoilz-1,0.0);
    std::vector<double> z_mid(nsoilz-1,0.0);
    std::vector<double> r(nsoilz-1,0.0);
    std::vector<double> e(nsoilz-1,0.0);
    std::vector<double> f(nsoilz-1,0.0);
    std::vector<double> g(nsoilz-1,0.0);
    std::vector<double> scalar(nsoilz);
    
    // compute soil or moisture conductivity
    if (type==1) {
        std::cout<<"Solving heat diffusion"<<std::endl;
        K_n = soil::soilThermalTransfer(psi_nsat,porosity,soil_q,b,Ci,nsoilz,1);
    }
    if (type==2) {
        std::cout<<"Solving moisture diffusion"<<std::endl;
        std::tie(D_n,K_n) = soil::soilMoistureTransfer(psi_nsat,K_nsat,porosity,soil_q,b,nsoilz);    
    }
    
    // here we set up and solve a tridiagonal matrix
    // AT(n+1) = r(n), where n denotes the time level
    // e, f, g the components of A matrix
    // T(n+1)  the soil temperature vector at t=n+1
    // r(n)    the soil temperature vector at t=n multiplied by coefficients         
    
    // interpolate soil_z and K_n to mid-points between 
    if (type==1) {
        for (int i=0; i<nsoilz-1; i++) {
            K_mid[i] = 0.5*(K_n[i]+K_n[i+1]);
            z_mid[i] = 0.5*(soil_z[i]+soil_z[i+1]);
        }
        scalar = soil_T;
        surf_scalar_last = surf_T_last;
        
    } else {
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
            
        if (type==2) {
            
            dz_p = soil_z[i] - soil_z[i+1];
            dz_m = soil_z[i+1] - soil_z[i+2];
            a_m  = -dz_p / (dz_m * (dz_p + dz_m));
            a_o  = (dz_p - dz_m) / (dz_p * dz_m);
            a_p  = dz_m / (dz_p * (dz_p + dz_m) );
            
            dKdz = dt*(K_n[i]*a_p + K_n[i+1]*a_o + K_n[i+2]*a_m);
                     
//            dKdz = dt*(K_n[i]*(soil_z[i+1]-soil_z[i+2])
//                   / ((soil_z[i]-soil_z[i+1])*(soil_z[i]-soil_z[i+2]))
//                 + K_n[i+1]*(2*soil_z[i+1]-soil_z[i]-soil_z[i+2])
//                   / ((soil_z[i+1]-soil_z[i])*(soil_z[i+1]-soil_z[i+2]))
//                 + K_n[i+2]*(soil_z[i+1]-soil_z[i])
//                   / ((soil_z[i+2]-soil_z[i])*(soil_z[i+2]-soil_z[i+1])));
            
            r[i] = r[i] + dKdz;
        }
    }
    
    // matrix coefficients for bottom level
    // first, construct ghost values
    int j       = nsoilz-2;
    double z_g  = 2*soil_z[j+1] - soil_z[j];
    double z_mg = (soil_z[j+1] + z_g) / 2.;
    
    C_p  = dt*K_mid[j]/(2*(z_mid[j]-z_mg)*(soil_z[j+1]-soil_z[j]));
    C_m  = dt*K_mid[j]/(2*(z_mid[j]-z_mg)*(soil_z[j]-z_g));
    C_tp = (1 + C_p + C_m);
    C_tm = (1 - C_p - C_m);
    
    e[j] = C_m  - C_p;
    f[j] = C_tp - C_m;
    r[j] = (C_p - C_m)*scalar[j] + (C_tm+C_m)*scalar[j+1];
    
    //scalar[nsoilz-2]-(scalar[nsoilz-2]-scalar[nsoilz-1])*dt*K_mid[nsoilz-2]
    //                               /(2*std::pow(soil_z[nsoilz-2]-soil_z[nsoilz-1],2));
    if (type==2) {
        dKdz = dt*(K_n[nsoilz-2]-K_n[nsoilz-1])/(soil_z[nsoilz-2]-soil_z[nsoilz-1]);
        r[nsoilz-2] = r[nsoilz-2] + dKdz;
    }
    
    // solve the tridiagonal system
    for (int kk=0; kk<nsoilz; kk++)
        std::cout<<"Before "<<scalar[kk]<<std::endl;
    try {
        if (type==1) matrix::tridiagonal(e,f,g,r,soil_T);
        if (type==2) matrix::tridiagonal(e,f,g,r,soil_q);
    } catch(std::string &e) {
        std::cout<<e<<std::endl;
        std::exit(0);
    }
    for (int kk=0; kk<nsoilz; kk++)
        std::cout<<"After "<<scalar[kk]<<std::endl;
}