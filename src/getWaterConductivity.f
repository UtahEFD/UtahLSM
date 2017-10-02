      subroutine getWaterConductivity(moisture,diffCond,hydrCond,
     +     porosity,satPotential,satHydrCond,soilExponent)
      
      implicit none
      
      real*8, dimension(:):: moisture,diffCond,hydrCond,porosity,
     +     satPotential,satHydrCond,soilExponent

!     Emperical relationships from Clapp and Hornberger (1978)                                                                                                 
!     need to be non-dimensionalized                                

      diffCond = -(soilExponent(:)*satHydrCond(:)*
     >     satPotential(:)/moisture(:))
     >     *(moisture(:)/porosity(:))**(soilExponent(:)+3.d0)
      
      hydrCond  = satHydrCond(:)*(moisture(:)/porosity(:))**
     >     (2.d0*soilExponent(:) + 3.d0)

      end subroutine getWaterConductivity
