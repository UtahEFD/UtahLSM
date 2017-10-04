      subroutine getWaterConductivity(moisture,diffCond,hydrCond,
     +     porosity1,satPotential1,satHydrCond1,soilExponent1)
      
      implicit none
      
      real*8, dimension(:):: moisture,diffCond,hydrCond,porosity1,
     +     satPotential1,satHydrCond1,soilExponent1

!     Emperical relationships from Clapp and Hornberger (1978)                                                                                                 
!     need to be non-dimensionalized                                

      diffCond = -(soilExponent1(:)*satHydrCond1(:)*
     >     satPotential1(:)/moisture(:))
     >     *(moisture(:)/porosity1(:))**(soilExponent1(:)+3.d0)
     
      hydrCond  = satHydrCond1(:)*(moisture(:)/porosity1(:))**
     >     (2.d0*soilExponent1(:) + 3.d0)

      end subroutine getWaterConductivity