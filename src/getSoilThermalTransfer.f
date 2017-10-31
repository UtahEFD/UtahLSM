      subroutine getSoilThermalTransfer(flag,moisture,porosity1,
     +     satPotential1,soilExponent1,heatCapSoil1,thermalTransfer)
      
      use globals
      use SEBmodule
      implicit none
      
      integer*4 flag
      real*8, dimension(:):: moisture,thermalTransfer,porosity1,
     +     satPotential1,soilExponent1,heatCapSoil1
      
      real*8, dimension(size(moisture)):: Pf, heatCap, moistPotential
      integer i
      
      ! flag = 0, thermal conductivity [Q = v dT/dz]                  
      ! flag = 1, thermal diffusivity [dT/dt = d( kdT/dz)dz]                

      ! compute soil conductivity from emperical formula (McCumber 1980)                 
      ! sat pot is expected in [cm] for Pf, so we multiply by 100
      ! resulting soil conductivity is in [ J/(m*s*K) ]
      moistPotential = 100.0d0*satPotential(:)*
     >     (porosity(:)/moisture(:))**soilExponent(:)
      
      Pf = log10( abs(moistPotential) )
      
      do i=1,size(Pf)
         if ( Pf(i).le.(5.1) )then
            thermalTransfer(i) = 418.46d0*exp( -( Pf(i)+2.7d0 ) )
         else
            thermalTransfer(i) = 0.172d0
         endif
      enddo
      
      ! convert to thermal diffusivity if flag=1   
      if(flag==1)then
         heatCap = (1-porosity1(:))*heatCapSoil1(:)
     >        + moisture*heatCapWater
         thermalTransfer = thermalTransfer/heatCap
      endif
      
      end subroutine getSoilThermalTransfer