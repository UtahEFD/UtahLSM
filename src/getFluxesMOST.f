      subroutine getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
     +     psi,psi0,fi,psiH,psiH0,fiH)     
      
      use globals
      use SEBmodule
      implicit none
      
      real*8  ustar,Uref,psi,psi0,fi,psiH,psiH0,fiH
      real*8, dimension(:) :: scalarFlux, scalarRef
      
      integer*4 i
      real*8 atmTemp,tempFlux,mixingRatio,obukhovLength,
     +     x,y,denom,q_gnd
      
      ! atmospheric scalars
      atmTemp     = scalarRef(temperatureIndex)
      mixingRatio = scalarRef(moistureIndex)
      
      ! iterate to solve for u* and scalar fluxes
      do i=1,4
                  
         ! compute friction velocity
         denom = dlog( z_m / z_o ) - psi + psi0
         ustar = Uref*vonk/denom
                  
         ! compute heat flux
         denom = dlog( z_s / z_t ) - psiH + psiH0 
         scalarFlux(temperatureIndex) = 
     >        ( gndScalars( 1,temperatureIndex )
     >        - scalarRef(temperatureIndex) ) * ustar*vonk/denom
         
         ! compute latent heat flux           
         call getSurfaceMixingRatio(q_gnd)
         
         scalarFlux(moistureIndex) = ( q_gnd 
     >   - scalarRef(moistureIndex) )*ustar*vonk/denom

         ! compute psi and fi functions for momentum and scalars 
         tempFlux = atmTemp*0.61d0*scalarFlux(moistureIndex)
     >            + scalarFlux(temperatureIndex)
     >            * (1.d0+0.61d0*mixingRatio)
         
         obukhovLength = -ustar**3*atmTemp/(vonk*grav*tempFlux)
         
         print*, "z/L=  ", z_m/obukhovLength
         print*, "u*= ", ustar
         
        ! the cap can probably be lifted here when not LES !!!                                                                                                               
!         if ( z_m/obukhovLength .gt. 5.0 )then
!            obukhovLength = z_m/5.0
!         elseif( z_m/obukhovLength .lt. -5.0)then
!            obukhovLength = -z_m/5.0
!         endif
        
         ! unstable 
         if( tempFlux.gt.0.0 )then
            
            ! momentum         
            x=(1.d0-(16.d0*z_m/obukhovLength))**0.25d0
            y=(1.d0-(16.d0*z_o/obukhovLength))**0.25d0           
            psi = 2.d0*dlog(0.5d0*(1+x)) +
     >            dlog(0.5d0*(1.d0+x**2.d0))-2.d0*atan(x)+pi/2.d0
            psi0 = 2.d0*dlog(0.5d0*(1+y)) -
     >            dlog(0.5d0*(1.d0+y**2.d0))-2.d0*atan(y)+pi/2.d0
            fi = 1.d0/x
	    
            ! scalar    
            x=(1.d0-(16.d0*z_s/obukhovLength))**0.25d0
            y=(1.d0-(16.d0*z_t/obukhovLength))**0.25d0                                                                                                                                              
            psiH  = 2.d0*dlog(0.5d0*(1.d0+x**2.d0))
            psiH0 = 2.d0*dlog(0.5d0*(1.d0+y**2.d0))
            fiH   = fi**2.d0
	    
         ! stable
         elseif( tempFlux.lt.-0.0 )then
            
            ! momentum                                                                                                                                     
            psi   = -5.d0*z_m/obukhovLength
            psi0  = -5.d0*z_o/obukhovLength                                                                                                                
            fi    = 1.d0 - psi
            
            ! scalar
            psiH  = -5.d0*z_s/obukhovLength
            psiH0 = -5.d0*z_t/obukhovLength                                    
            fiH   = 0.74d0 - psiH
         
         ! neutral
         else
            
            ! momentum
            psi   = 0.d0
            psi0  = 0.d0
            fi    = 1.d0
            
            ! scalar
            psiH  = psi
            psiH0 = psi0
            fiH   = 0.74d0
            
         endif 
            
      enddo
      
      end subroutine getFluxesMOST