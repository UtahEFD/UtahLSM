      subroutine getWaterConductivity(moisture,diffCond,hydrCond,
     +     por,satPot,satHCond,soilEx)
      
      implicit none
      
      real*8, dimension(:):: moisture,diffCond,hydrCond,por,
     +     satPot,satHCond,soilEx

!     Emperical relationships from Clapp and Hornberger (1978)                                                                                                 
!     need to be non-dimensionalized                                

      diffCond = -(soilEx(:)*satHCond(:)*satPot(:)/moisture(:))
     >     *(moisture(:)/por(:))**(soilEx(:)+3.d0)
     
      hydrCond  = satHCond(:)*(moisture(:)/por(:))**
     >     (2.d0*soilEx(:) + 3.d0)

      end subroutine getWaterConductivity