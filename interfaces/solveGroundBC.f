      subroutine solveGroundBC(Uref,scalarRef, 
     +     gndScalars,ustar,scalarFlux,soilHeatFlux,porosity,
     +     satPotential,satHydrCond,soilExponent,heatCapSoil,albedo,
     +     minAlbedo,measRad,netRad,Psi,Psi0,fi,fiH,zlevel,zo,zGnd)
      
      real*8 Uref,ustar,soilHeatFlux,albedo,minAlbedo,
     +     netRad,zo,Psi,Psi0,fi,fiH,zlevel
      real*8,dimension(:) :: scalarRef,scalarFlux,measRad,
     +     zGnd,porosity,satPotential,satHydrCond,soilExponent,
     +     heatCapSoil
      real*8,dimension(:,:) :: gndScalars
      end subroutine solveGroundBC
