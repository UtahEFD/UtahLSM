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
                 std::vector<double>& soil_T,std::vector<double>& soil_T_last,
                 std::vector<double>& soil_q,std::vector<double>& soil_q_last,
                 int julian_day, double utc, double latitude, double longitude,
                 double albedo, double emissivity, double R_net, int comp_rad,
                 double& zeta_m,double& zeta_s,double& zeta_o,double& zeta_t,
                 double& ustar, double& flux_wT, double& flux_wq, double& flux_gr) : 
                 first(first),dt(dt),z_o(z_o),z_t(z_t),z_m(z_m),z_s(z_s),
                 atm_p(atm_p),atm_ws(atm_ws),atm_T(atm_T),atm_q(atm_q),
                 nsoilz(nsoilz),soil_z(soil_z),soil_type(soil_type),soil_T(soil_T),
                 soil_T_last(soil_T_last),soil_q(soil_q),soil_q_last(soil_q_last),
                 julian_day(julian_day),utc(utc),latitude(latitude),longitude(longitude), 
                 albedo(albedo),emissivity(emissivity),R_net(R_net),comp_rad(comp_rad),
                 zeta_m(zeta_m),zeta_s(zeta_s),zeta_o(zeta_o),zeta_t(zeta_t),
                 ustar(ustar),flux_wT(flux_wT),flux_wq(flux_wq),flux_gr(flux_gr) {
    
    // save incoming temp and moisture
    surf_T_last = soil_T[0];
    surf_q_last = soil_q[0];
    
    // set soil properties at each depth
    setSoilProperties();
    
    // initial time requires estimate of fluxes
    if (first) {
	    computeFluxes(soil_T[0],soil_q[0]);
	}
    
    // check whether users provide net radiation
    if (comp_rad==true) computeRadiation();
    
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
}

// Set soil properties at each depth
void UtahLSM :: setSoilProperties() {
	
	soil::soilProperties soilProperties = soil::soilTypeProperties(soil_type,nsoilz);

	b        = soilProperties.b;
	psi_sat  = soilProperties.psi_sat;
	porosity = soilProperties.porosity;
	K_sat    = soilProperties.K_sat;
	Ci       = soilProperties.Ci;
}

// compute surface fluxes using Monin-Obukhov w/Dyer-Hicks
void UtahLSM :: computeFluxes(double sfc_T, double sfc_q) {
    
    // local variables
    int max_iterations = 200;
    bool converged = false;
    double L=0, gnd_q, zot_q, zot_T, flux_wTv, atm_Tv;
    double last_L, criteria = 0.1, ref_T = 300.;
    
    // compute surface mixing ratio
    gnd_q  = soil::surfaceMixingRatio(psi_sat[0],porosity[0],b[0],sfc_T,sfc_q,atm_p);    
    atm_Tv = atm_T * (1 + 0.61*atm_q);
    
    // iterate for convergence
    for (int i=0; i<max_iterations; ++i) {
        
        // extrapolate ground temperature/moisture to roughness height
        zot_T = sfc_T + (atm_T-sfc_T)*std::log(z_o/z_t)*most::fh(z_s/z_t,zeta_s,zeta_t)/c::vonk;
        zot_q = gnd_q + (atm_q-gnd_q)*std::log(z_o/z_t)*most::fh(z_s/z_t,zeta_s,zeta_t)/c::vonk;
         
        //zot_T = sfc_T;
        //zot_q = gnd_q;
        
        // compute friction velocity
        ustar = atm_ws*most::fm(z_m/z_o,zeta_m,zeta_o);
        
        // compute heat flux
        flux_wT = (zot_T-atm_T)*ustar*most::fh(z_s/z_o,zeta_s,zeta_o);
        
        // compute latent heat flux
	    //if (first){
		//    flux_wq   = (c::Cp_air/c::Lv)*flux_wT/(-4);
		//    soil_q[0] = soil::surfaceSoilMoisture(psi_sat[0],porosity[0],b[0],
        //                                         sfc_T,gnd_q,atm_p);
		//} else {
		flux_wq = (zot_q-atm_q)*ustar*most::fh(z_s/z_o,zeta_s,zeta_o);
	    //} 
        
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
        if (converged) break;
    }
        
    if (!converged) {
	    std::cout<<"Something is wrong: L didn't converge"<<std::endl;
	    std::cout<<flux_wq<<" "<<soil_q[0]<<std::endl;
	    throw(1);
	}
}

// compute net radiation
void UtahLSM :: computeRadiation() {

    R_net = radiation::net(soil_T[0],emissivity,julian_day,utc,
                           latitude,longitude,albedo);
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
    double flux_criteria = 0.0001;
            
    // compute SEB using current bracketed temperatures
    SEB_l = computeSEB(temp_1);
    SEB_h = computeSEB(temp_2);
	
    // dynamic bracket adjustments
    bool out_of_bracket = (SEB_l > 0.0 && SEB_h > 0.0) || (SEB_l < 0.0 && SEB_h < 0.0);
	
    while (out_of_bracket) {

	    // expand brackets by 5 K
	    temp_1 -= 5;
	    temp_2 += 5;
	    
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
            if (std::abs(dTs) <= temp_criteria) break;
            
            // if convergence fails, recompute flux
            computeFluxes(soil_T[0],soil_q[0]);
        }
        
        // save current flux for convergence criteria
        last_F = flux_wT;
        
        // recompute heat flux using new temperature
        computeFluxes(soil_T[0],soil_q[0]);
        
        // check for convergence
        if (std::abs(flux_wT-last_F) <= flux_criteria) {
	        std::cout<<std::endl;
	        SEB = computeSEB(soil_T[0]);
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

    
    // compute heat storage term
    K_soil     = soil::soilThermalTransfer(psi_sat,porosity,soil_q,b,Ci,depth,0);

    if (first) {
    	//Qg =  0.5*(K_soil[0]+K_soil[1])*(sfc_T - soil_T[1])/(soil_z[0]-soil_z[1]);
    	//Qg =  K_soil[1]*(sfc_T - soil_T[2])/(soil_z[0]-soil_z[2]);
    	Qg = 0.4*R_net;
    } else {
        // compute heat flux within soil and find depth of minimum
        // this is the depth to integrate time change of T		
        int int_depth = 0;
        double hfs;
        double min_sflux = 10000.;
        for (int d=0; d<nsoilz-1; ++d) {
            
            hfs = (d==0) ? 0.5*(K_soil[d]+K_soil[d+1]) * (sfc_T     - soil_T[d+1])/(soil_z[d]-soil_z[d+1]):
                           0.5*(K_soil[d]+K_soil[d+1]) * (soil_T[d] - soil_T[d+1])/(soil_z[d]-soil_z[d+1]);
            if (std::abs(hfs)<std::abs(min_sflux)) {
                min_sflux = hfs;
                int_depth = d+1;
            }
        }
            	
        Qg = 0;
        //int_depth = 1;
        // Integrate time change to depth of minimum flux
        for (int d=0; d<int_depth; ++d) {			
            double heat_cap1 = (1-porosity[d  ])*Ci[d  ] + soil_q[d  ]*c::Ci_wat + (porosity[d  ]-soil_q[d  ])*c::Cp_air;
            double heat_cap2 = (1-porosity[d+1])*Ci[d+1] + soil_q[d+1]*c::Ci_wat + (porosity[d+1]-soil_q[d+1])*c::Cp_air;
            double dz  = soil_z[d]-soil_z[d+1];
            double dT1 = (d==0) ? sfc_T-soil_T_last[d] : soil_T[d]-soil_T_last[d];
            double dT2 = soil_T[d+1]-soil_T_last[d+1];
            
            Qg = Qg + 0.5*( heat_cap1*dT1 + heat_cap2*dT2)*(dz/dt);
               
        }
        
        //std::cout<<Qg<<std::endl;
        //Qg = Qg + 0.5*(K_soil[0]+K_soil[1])*(sfc_T - soil_T[1])/(soil_z[0]-soil_z[1]);
        //Qg = Qg + K_soil[1]*(sfc_T - soil_T[2])/(soil_z[0]-soil_z[2]);
        //std::cout<<Qg<<std::endl;
    }    
    flux_gr    = Qg;
    
    
    // compute surface energy balance
    SEB = R_net - Qg - Qh - Ql;    
    return SEB;
}

// computes derivative of SEB wrt surface temperature
double UtahLSM :: computeDSEB(double sfc_T) {
    
    // local variables
    double heat_cap, dSEB_dT;
    
    // compute heat capacity of soil surface
    heat_cap = (1-porosity[0])*Ci[0] + soil_q[0]*c::Ci_wat + (porosity[0]-soil_q[0])*c::Cp_air;
    
    //compute derivative of SEB wrt temperature
    dSEB_dT = - 4.*emissivity*c::sb*std::pow(sfc_T,3.) 
              - c::rho_air*c::Cp_air*ustar*most::fh(z_s/z_t,zeta_s,zeta_t)
              - heat_cap*(soil_z[0]-soil_z[1])/(2*dt);
    return dSEB_dT;
}


// Solve the surface energy balance
void UtahLSM :: solveSMB() {
        
    // local variables
    int max_iter_temp = 200;
    int max_iter_flux = 200; 
    double dqs, dqs_old;
    double SMB, dSMB_dq, SMB_l, SMB_h;
    double mois_l, mois_h, last_q, last_F;
    double mois_1 = soil_q[0] - 0.01;
    double mois_2 = soil_q[0] + 0.01;
    double mois_criteria = 0.001;
    double flux_criteria = 0.000001;
            
    // compute SMB using current bracketed surface moisture
    SMB_l = computeSMB(mois_1);
    SMB_h = computeSMB(mois_2);
	
    // dynamic bracket adjustments
    bool out_of_bracket = (SMB_l > 0.0 && SMB_h > 0.0) || (SMB_l < 0.0 && SMB_h < 0.0);
	
    while (out_of_bracket) {

	    // expand brackets by 0.05
	    mois_1 -= 0.05;
	    mois_2 += 0.05;
	    
	    // recompute SMB at brackets
	    SMB_l = computeSMB(mois_1);
		SMB_h = computeSMB(mois_2);
		
		// check for proper brackets 
		out_of_bracket = (SMB_l > 0.0 && SMB_h > 0.0) || (SMB_l < 0.0 && SMB_h < 0.0);
    }
    
    if ((SMB_l > 0.0 && SMB_h > 0.0) || (SMB_l < 0.0 && SMB_h < 0.0)) {
        throw("Please adjust brackets for qs");
    }
    
    // if SEB from low bracket qs = 0, then that value of qs is solution
    if (SMB_l == 0.0) soil_q[0] = mois_1;

    // if SEB from high bracket qs = 0, then that value of qs is solution
    if (SMB_h == 0.0) soil_q[0] = mois_2;
    
    // orient the solutions such that SEB(temp_l) < 0;
    if (SMB_l < 0.0) {
        mois_l = mois_1;
        mois_h = mois_2;
    } else {
        mois_l = mois_2;
        mois_h = mois_1; 
    }
        
    // prepare for convergence looping
    dqs = std::abs(mois_h-mois_l);
     
    // convergence loop for flux
    for (int ff = 0; ff < max_iter_flux; ff++) {
        
        // convergence loop for moisture
        for (int tt = 0; tt < max_iter_temp; tt++) {
            
            // compute SMB and dSMB_dqs
            SMB     = computeSMB(soil_q[0]);
            dSMB_dq = computeDSMB(soil_q[0]);
                   
            // update brackets
            if (SMB<0.) mois_l = soil_q[0];
            if (SMB>0.) mois_h = soil_q[0];
            
            // bracket and bisect temperature if Newton out of range
            if ((((soil_q[0]-mois_h)*dSMB_dq-SMB)*((soil_q[0]-mois_l)*dSMB_dq-SMB)>0.0) 
               || (std::abs(2*SMB) > std::abs(dqs_old*dSMB_dq))) {
                dqs_old   = dqs;
                dqs       = 0.5*(mois_h-mois_l);
                last_q    = soil_q[0];
                soil_q[0] = mois_l + dqs;
                if (mois_l == soil_q[0]) break;
            } else {
                dqs_old   = dqs;
                dqs       = SMB / dSMB_dq;
                last_q    = soil_q[0];
                soil_q[0] = soil_q[0] - dqs;
                if (last_q == soil_q[0]) break;
            }
            
            // check for convergence
            if (std::abs(dqs) <= mois_criteria) break;
            
            // if convergence fails, recompute flux
            computeFluxes(soil_T[0],soil_q[0]);
        }
        
        // save current flux for convergence criteria
        last_F = flux_wq;
        
        // recompute heat flux using new temperature
        computeFluxes(soil_T[0],soil_q[0]);
        
        // check for convergence
        if (std::abs(flux_wq-last_F) <= flux_criteria) {
	        SMB = computeSMB(soil_q[0]);
	        double H = c::rho_air*c::Cp_air*flux_wT;
            double L = c::rho_air*c::Lv*flux_wq;
            double G = flux_gr;
            std::cout<<"--- Converged with heat flux = "<<H<<", mois flux = "<<L<<", grnd flux = "<<G<<", and temp = "<<soil_T[0]<<std::endl;
	        break;
	    }
                
        // if flux fails to converge, split moisture
        soil_q[0] = 0.5*(soil_q[0] + last_q);
        
        // recompute stability corrections and fluxes
        computeFluxes(soil_T[0],soil_q[0]);
    }
}

// computes SMB
double UtahLSM :: computeSMB(double sfc_q) {
    
    // local variables
    int depth = nsoilz;
    double W, E, P, SMB;
    double psi_n0, psi_n1, psi_n2;
    soil::soilTransfer transfer;
    std::vector<double> K_n;
    std::vector<double> q_prof;
    
    // Set precipitation to zero for now
    P = 0;
    
    // adjust surface value for transfer calculation
    q_prof    = soil_q;
    q_prof[0] = sfc_q;
    
    // moisture potential at first level below ground
    psi_n0 = std::abs(psi_sat[0])*std::pow(porosity[0]/sfc_q,b[0]);
    psi_n1 = std::abs(psi_sat[1])*std::pow(porosity[1]/soil_q[1],b[1]);
    psi_n2 = std::abs(psi_sat[2])*std::pow(porosity[2]/soil_q[2],b[2]);
    
    // compute initial soil moisture flux
    transfer = soil::soilMoistureTransfer(psi_sat,K_sat,porosity,q_prof,b,depth);
    K_n      = transfer.k;
        
    if (first) {
    	//W = c::rho_wat*0.5*(K_n[0]+K_n[1])*((psi_n0 - psi_n1)/(soil_z[0]-soil_z[1]) + 1);
    	W = c::rho_wat*K_n[1]*((psi_n0 - psi_n2)/(soil_z[0]-soil_z[2]) + 1);
    } else {
        // compute moisture flux within soil and find depth of minimum
        // this is the depth to integrate time change of q		
        int int_depth = 0;
        double mfs;
        double min_mflux = 10000.;
        for (int d=0; d<nsoilz-1; ++d) {
            
            // moisture potential at first level below ground
			psi_n0 = (d==0) ? std::abs(psi_sat[d])*std::pow(porosity[d]/sfc_q,    b[d]):
			                  std::abs(psi_sat[d])*std::pow(porosity[d]/soil_q[d],b[d]);
			psi_n1 = std::abs(psi_sat[d+1])*std::pow(porosity[d+1]/soil_q[d+1],b[d+1]);
            mfs    = 0.5*c::rho_wat*(K_n[d]+K_n[d+1])*((psi_n0 - psi_n1)/(soil_z[0]-soil_z[1]) + 1);
            
            if (std::abs(mfs)<std::abs(min_mflux)) {
                min_mflux = mfs;
                int_depth = d+1;
            }
        }
            	
        W = 0;
        //int_depth=1;
        // Integrate time change to depth of minimum flux
        for (int d=0; d<int_depth; ++d) {			
            
            double dz  = soil_z[d]-soil_z[d+1];
            double dq1 = (d==0) ? sfc_q-soil_q_last[d] : soil_q[d]-soil_q_last[d];
            double dq2 = soil_q[d+1]-soil_q_last[d+1];
            
            W = W + 0.5*c::rho_wat*(dq1 + dq2)*(dz/dt);
        }
        
        psi_n0 = std::abs(psi_sat[0])*std::pow(porosity[0]/sfc_q,b[0]);
		psi_n1 = std::abs(psi_sat[1])*std::pow(porosity[1]/soil_q[1],b[1]);
		psi_n2 = std::abs(psi_sat[2])*std::pow(porosity[2]/soil_q[2],b[2]);
        //W = W + c::rho_wat*0.5*(K_n[0]+K_n[1])*((psi_n0 - psi_n1)/(soil_z[0]-soil_z[1]) + 1);
        W = W + c::rho_wat*K_n[1]*((psi_n0 - psi_n2)/(soil_z[0]-soil_z[2]) + 1);
    }
    
    // update evaporation
    computeFluxes(soil_T[0],sfc_q);
    E = c::rho_air*flux_wq;
    //gnd_q = soil::surfaceMixingRatio(psi_sat[0],porosity[0],b[0],soil_T[0],sfc_q,atm_p); 
    //E     = c::rho_air*(gnd_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
    
    // compute surface energy balance
    SMB = W + E - P;    
    return SMB;
}

// computes derivative of SMB wrt surface moisture
double UtahLSM :: computeDSMB(double sfc_q) {
	
    
    // local variables
    double A, C, es, qs, dSMB_dT;
    
    es = 6.1078*std::exp(17.269*(soil_T[0]-273.15)/(soil_T[0]-35.86));
    qs = 0.622*(es/(atm_p-0.378*es));
	A  = c::rho_air*ustar*qs*most::fh(z_s/z_t,zeta_s,zeta_t);
	C  = (c::grav*psi_sat[0]*std::pow(porosity[0],b[0])) / (c::Rv*soil_T[0]);
    
    //compute derivative of SEB wrt temperature
    dSMB_dT = c::rho_wat*(soil_z[0]-soil_z[1])/(2*dt)
            - A*b[0]*C*std::pow(sfc_q,(-b[0]-1))*std::exp(C*sfc_q-b[0]);
    return dSMB_dT;
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
    
    // set up and solve a tridiagonal matrix
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
        D_n      = transfer.d;
        K_n      = transfer.k;
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
