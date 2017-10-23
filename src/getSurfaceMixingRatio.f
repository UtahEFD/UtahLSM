      subroutine getSurfaceMixingRatio(q_gnd)
      ! solves for surface mixiing ratio following McCumber (1980)
      
      use globals
      use SEBmodule
      implicit none
      
      real*8 q_gnd
      real*8 moistPotential(2), h, partialPressure, satHum, specHum_gnd
      
      ! McCumber (1980); see Pielke (2013), page 353                                                             
      moistPotential(1) = satPotential(1)*
     >     (porosity(1)/gndScalars(1,moistureIndex))**soilExponent(1)
      
      h = exp(g_hat*moistPotential(1)
     >     / (waterGasConst*gndScalars(1,temperatureIndex)))
      
      partialPressure = 6.1078d0*exp(17.269d0*
     >     (gndScalars(1,temperatureIndex) - 273.16d0/scalarScales(1))
     >     / (gndScalars(1,temperatureIndex) - 35.86d0/scalarScales(1)))
      
      satHum = 0.622d0*( partialPressure /
     >     ( pressureScale - 0.378d0*partialPressure) )
      
      specHum_gnd = h*satHum/scalarScales(2)
      
      ! convert from specific humidity to mixing ratio                                                                                                                         
      q_gnd = specHum_gnd/(1-specHum_gnd)
      
      end subroutine getSurfaceMixingRatio