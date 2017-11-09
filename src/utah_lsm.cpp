//
//  utah_lsm.cpp
//  
//  This class handles the Utah LSM
//
//  Created by Jeremy Gibbs on 10/31/17.
//

#include "utah_lsm.hpp"
#include "constants.hpp"
#include "soil.hpp"
#include "most.hpp"
#include <cmath>
#include <iostream>

namespace {
    namespace consts = Constants;
}

UtahLSM::UtahLSM(double z_o,double z_t,double z_m,double z_s, 
                double atm_p,double atm_ws,double atm_T,double atm_q, 
                int nsoilz,double* soil_z,double* soil_T,double* soil_q,
                double* porosity,double* psi_nsat,double* K_nsat,double* b,double* Ci,
                int julian_day, double utc, double latitude, double longitude,
                double albedo, double emissivity, double R_net,
                double& zeta_m,double& zeta_s,double& zeta_o,double& zeta_t,
                double& ustar, double& flux_wT, double& flux_wq) : 
                z_o(z_o),z_t(z_t),z_m(z_m),z_s(z_s),
                atm_p(atm_p),atm_ws(atm_ws),atm_T(atm_T),atm_q(atm_q),
                nsoilz(nsoilz),julian_day(nsoilz),utc(utc),latitude(latitude), 
                longitude(longitude), albedo(albedo),emissivity(emissivity),R_net(R_net),
                soil_z(soil_z),soil_T(soil_T),soil_q(soil_q),
                porosity(porosity),psi_nsat(psi_nsat), K_nsat(K_nsat),b(b),Ci(Ci),
                zeta_m(zeta_m),zeta_s(zeta_s),zeta_o(zeta_o),zeta_t(zeta_t),
                ustar(ustar),flux_wT(flux_wT),flux_wq(flux_wq) {
    
    // run the model    
    computeFluxes();
    solveSEB();
}

// compute surface fluxes using Monin-Obukhov w/Dyer-Hicks
void UtahLSM :: computeFluxes() {
    
    // local variables
    double L, sfc_q, flux_wTv;
    
    // compute surface mixing ratio
    sfc_q = soil::surfaceMixingRatio(psi_nsat[0],porosity[0],b[0],
                                     soil_T[0],soil_q[0],atm_p);
    
    // iterate to solve for u* and scalar fluxes
    for (int i=0; i<4; ++i) {
        
        // compute friction velocity
        ustar = atm_ws * most::fm(z_m/z_o,zeta_m,zeta_o);
        
        // compute heat flux
        flux_wT = (soil_T[0]-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);

        // compute latent heat flux
        flux_wq = (sfc_q-atm_q)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
        
        // compute virtual heat flux
        flux_wTv = flux_wT + atm_T*0.61*flux_wq;
        
        // compute L
        L = std::pow(-ustar,3.)*atm_T/(consts::vonk*consts::grav*flux_wTv);
        
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

void UtahLSM :: solveSEB() {
    
    fluxConverged = false;
    tempConverged = false;
    
    // save current surface temperature for initial bisection
    lastTemp = soil_T[0];
        
    // outer flux loop
    do {
        
        // inner temperature loop
        do {
            tempConverged = convergeTemp();
            
        } while (!tempConverged);        
        
        std::cout<<"Temp: "<<soil_T[0]<<" Flux: "<<flux_wT<<std::endl;
        fluxConverged = convergeFlux();
                         
    } while (!fluxConverged);
    
}

// use Newton-Raphson for temperature convergence
bool UtahLSM :: convergeTemp() {
    
    // local variables
    double Ksoil, Qg, Qh, Ql, SEB;
    double dQh_dT, dSEB_dT, dTs;
    double tempCriteria = 0.01; 
    
    // compute soil heat flux
    Ksoil = soil::soilThermalTransfer(psi_nsat,porosity,soil_q,b,Ci,2,0);
    Qg    = -(soil_T[1]-soil_T[0])*(Ksoil/(soil_z[1]-soil_z[0]));
    
    // write sensible and latent heat fluxes in [W/m^2]
    Qh = Constants::rho_air*consts::Cp_air*flux_wT;
    Ql = Constants::rho_air*consts::Lv*flux_wq;
    
    // compute surface energy balance
    SEB = Qg + Qh + Ql - R_net;
    
    // compute derivative of heat flux wrt temperature
    dQh_dT = ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
    dQh_dT = consts::rho_air*consts::Cp_air*dQh_dT;
    
    //compute derivative of SEB wrt temperature
    dSEB_dT = 4.*emissivity*consts::sb*std::pow(soil_T[0],3.)
              + Ksoil/(soil_z[1]-soil_z[0]) + dQh_dT;
    
    // compute time change in temperature
    dTs = SEB / dSEB_dT;
    //std::cout<<dTs<<std::endl;
    // update surface temperature
    soil_T[0] = soil_T[0] - dTs;
    
    // check for convergence
    return std::abs(dTs) <= tempCriteria;
}

// use converged temperature to acheive flux convergence
bool UtahLSM :: convergeFlux() {
    
    // local variables
    bool converged = false;
    double fluxCriteria = .001;
    
    // save current flux for convergence criteria
    lastFlux = flux_wT;
    
    // recompute heat flux using new temperature
    //flux_wT = (soil_T[0]-atm_T)*ustar*most::fh(z_s/z_t,zeta_s,zeta_t);
    computeFluxes();
    
    std::cout<<std::abs(flux_wT-lastFlux)<<std::endl;
    // check for convergence
    converged = std::abs(flux_wT-lastFlux) <= fluxCriteria;
    
    if (!converged) {
        // if not converged, bisect temperature
        soil_T[0] = (soil_T[0]+lastTemp)/2.;
        
        // save bisected temperature for next loop
        lastTemp = soil_T[0];
        
        // recompute fluxes with bisected temperature
        computeFluxes();
    }
    
    return converged;
}