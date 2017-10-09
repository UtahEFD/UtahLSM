      subroutine getSoilThermalTransfer(moisture,thermalTransfer,flag,
     +     porosity1,satPotential1,soilExponent1,heatCapSoil1)
      
      use globals
      use SEBmodule
      implicit none
      
      integer*4 flag
      real*8, dimension(:):: moisture,thermalTransfer,porosity1,
     +     satPotential1,soilExponent1,heatCapSoil1
      
      real*8, dimension(size(moisture)):: Pf, heatCap
      integer i
!     flag = 0 if thermal conductivity is requested (for solving Q = v dT/dz )                  
!     flag = 1 if thermal diffusivity is requested (for solving dT/dt = d( k dT/dz)dz                

!     compute soil conductivity from emperical formula (McCumber 1980)                 
      Pf = log10( abs(100.d0*z_i*satPotential1(:)*
     >     (porosity1(:)/moisture)**soilExponent1(:) ) )
      
      do i=1,size(Pf)
         if ( Pf(i).le.(5.1) )then
            thermalTransfer(i) = 418.46d0*exp( -( Pf(i)+2.7d0 ) )
         else
            thermalTransfer(i) = 0.172d0
         endif
      enddo
     
!     emperical relationship for thermal conductivity is in [ J/(m*s*K) ]   
!     non-dimensionalize the conductivity 
                                    
      if(flag==1)then
         heatCap = (1-porosity1(:))*heatCapSoil1(:)
     >        + moisture*heatCapWater
         thermalTransfer = thermalTransfer
     >        /(heatCap*densityAir*Cp_air*z_i*uScale)
      else
         thermalTransfer = thermalTransfer
     >        /(densityAir*Cp_air*z_i*uScale)
      endif
      
      end subroutine getSoilThermalTransfer

