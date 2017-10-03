      subroutine solveGroundBC(Uref,scalarRef, 
     +     gndScalars,ustar,scalarFlux,soilHeatFlux,porosity,
     +     satPotential,satHydrCond,soilExponent,heatCapSoil,
     +     measRad,netRad,Psi,Psi0,fi,fiH,zGnd)
      
      real*8 Uref,ustar,soilHeatFlux,netRad,Psi,Psi0,fi,fiH
      real*8,dimension(:) :: scalarRef,scalarFlux,measRad,
     +     zGnd,porosity,satPotential,satHydrCond,soilExponent,
     +     heatCapSoil
      real*8,dimension(:,:) :: gndScalars
      end subroutine solveGroundBC
