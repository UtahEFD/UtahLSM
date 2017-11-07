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
//#include <netcdf>
//using namespace netCDF;
//using namespace netCDF::exceptions;

UtahLSM::UtahLSM(double z_o,double z_t,double z_m,double z_s, 
                double atm_p,double atm_ws,double atm_T,double atm_q, 
                int nsoilz,double* soil_z,double* soil_T,double* soil_q,
                double* porosity,double* psi_nsat,double* K_nsat,double* b,double* Ci,
                int julian_day, double utc, double latitude, double longitude,
                double albedo, double emissivity, double R_net,
                double* phiM, double* psiM, double* psiM0,
                double* phiH, double* psiH, double* psiH0,
                double* ustar, double* flux_wT, double* flux_wq) : 
                z_o(z_o),z_t(z_t),z_m(z_m),z_s(z_s),
                atm_p(atm_p),atm_ws(atm_ws),atm_T(atm_T),atm_q(atm_q),
                nsoilz(nsoilz),julian_day(nsoilz),utc(utc),latitude(latitude), 
                longitude(longitude), albedo(albedo),emissivity(emissivity),R_net(R_net),
                soil_z(soil_z),soil_T(soil_T),soil_q(soil_q),
                porosity(porosity),psi_nsat(psi_nsat), K_nsat(K_nsat),b(b),Ci(Ci),
                phiM(phiM),psiM(psiM),psiM0(psiM0),phiH(phiH), psiH(psiH),psiH0(psiH0),
                ustar(ustar),flux_wT(flux_wT),flux_wq(flux_wq) {
    
    computeFluxes();
}

void UtahLSM :: computeFluxes() {
    
    double sfc_q, L, flux_wTv;
    
    // iterate to solve for u* and scalar fluxes
    for (int i=0; i<4; ++i) {
        
        // compute friction velocity
        *ustar   = atm_ws * Constants::vonk / (std::log(z_m/z_o) - (*psiM) + (*psiM0) );
        
        // compute heat flux
        *flux_wT = (soil_T[0]-atm_T)* (*ustar)*Constants::vonk/(std::log(z_m/z_o)- (*psiH)+(*psiH0));

        // compute latent heat flux
        sfc_q = soil::surfaceMixingRatio(psi_nsat[0],porosity[0],b[0],
                                         soil_T[0],soil_q[0],atm_p);
        *flux_wq = (sfc_q-atm_q)*(*ustar)*Constants::vonk/(std::log(z_s/z_t)- (*psiH) + (*psiH0));
        
        // compute virtual heat flux
        flux_wTv = atm_T*0.61*(*flux_wq) + (*flux_wT)*(1+0.61*sfc_q);
        
        // compute L
        L = std::pow(-(*ustar),3.) * atm_T/(Constants::vonk*Constants::grav*flux_wTv);
        
        // compute stability correction functions
        *psiM  = most::psim(z_m/L);
        *psiM0 = most::psim(z_o/L);
        *psiH  = most::psih(z_s/L);
        *psiH0 = most::psih(z_t/L);
        
    }
}
        
        
        
        // compute heat flux
        //*flux_wT = (soil_T[0]-atm_T)*ustar*Constants::vonk/(std::log(z_m/z_o)-psiH+psiH0)
        
        //compute
        
        //*flux_wq = (soil_q[0]-atm_q)*ustar*Constants::vonk/(std::log(z_s/z_t)-psiH+psiH0)
        
        // get surface mixing ratio
        
        //flux_wTv = atm_T*0.61*
            
    // iterate to solve for u* and scalar fluxes
      //for (int i=0; i<4; ++i) {
                  
         // compute friction velocity
        // * most::fmn(z_m,z_o);
                  
//         ! compute heat flux
//         denom = dlog( z_s / z_t ) - psiH + psiH0 
//         scalarFlux(temperatureIndex) = 
//     >        ( gndScalars( 1,temperatureIndex )
//     >        - scalarRef(temperatureIndex) ) * ustar*vonk/denom
//         
//         ! compute latent heat flux           
//         call getSurfaceMixingRatio(q_gnd)
//         
//         scalarFlux(moistureIndex) = ( q_gnd 
//     >   - scalarRef(moistureIndex) )*ustar*vonk/denom
//
//         ! compute psi and fi functions for momentum and scalars 
//         tempFlux = atmTemp*0.61d0*scalarFlux(moistureIndex)
//     >            + scalarFlux(temperatureIndex)
//     >            * (1.d0+0.61d0*mixingRatio)
//         
//         obukhovLength = -ustar**3*atmTemp/(vonk*grav*tempFlux)
//         
//         print*, "z/L=  ", z_m/obukhovLength
//         print*, "u*= ", ustar
//         
//        ! the cap can probably be lifted here when not LES !!!                                                                                                               
//!         if ( z_m/obukhovLength .gt. 5.0 )then
//!            obukhovLength = z_m/5.0
//!         elseif( z_m/obukhovLength .lt. -5.0)then
//!            obukhovLength = -z_m/5.0
//!         endif
//        
//         ! unstable 
//         if( tempFlux.gt.0.0 )then
//            
//            ! momentum         
//            x=(1.d0-(16.d0*z_m/obukhovLength))**0.25d0
//            y=(1.d0-(16.d0*z_o/obukhovLength))**0.25d0           
//            psi = 2.d0*dlog(0.5d0*(1+x)) +
//     >            dlog(0.5d0*(1.d0+x**2.d0))-2.d0*atan(x)+pi/2.d0
//            psi0 = 2.d0*dlog(0.5d0*(1+y)) -
//     >            dlog(0.5d0*(1.d0+y**2.d0))-2.d0*atan(y)+pi/2.d0
//            fi = 1.d0/x
//	    
//            ! scalar    
//            x=(1.d0-(16.d0*z_s/obukhovLength))**0.25d0
//            y=(1.d0-(16.d0*z_t/obukhovLength))**0.25d0                                                                                                                                              
//            psiH  = 2.d0*dlog(0.5d0*(1.d0+x**2.d0))
//            psiH0 = 2.d0*dlog(0.5d0*(1.d0+y**2.d0))
//            fiH   = fi**2.d0
//	    
//         ! stable
//         elseif( tempFlux.lt.-0.0 )then
//            
//            ! momentum                                                                                                                                     
//            psi   = -5.d0*z_m/obukhovLength
//            psi0  = -5.d0*z_o/obukhovLength                                                                                                                
//            fi    = 1.d0 - psi
//            
//            ! scalar
//            psiH  = -5.d0*z_s/obukhovLength
//            psiH0 = -5.d0*z_t/obukhovLength                                    
//            fiH   = 0.74d0 - psiH
//         
//         ! neutral
//         else
//            
//            ! momentum
//            psi   = 0.d0
//            psi0  = 0.d0
//            fi    = 1.d0
//            
//            ! scalar
//            psiH  = psi
//            psiH0 = psi0
//            fiH   = 0.74d0
//            
//         endif 
//            
//      enddo
   // }
//}