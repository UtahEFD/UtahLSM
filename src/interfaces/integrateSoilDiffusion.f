      subroutine integrateSoilDiffusion(gndScalars,lastSurfScalars,ind,
     +     zGnd,porosity,satPotential,soilExponent,heatCapSoil,
     +     satHydrCond)
      integer*4 ind
      real*8,dimension(:) :: lastSurfScalars,zGnd,porosity,satPotential,
     +     soilExponent,heatCapSoil,satHydrCond
      real*8,dimension(:,:) :: gndScalars
      end subroutine integrateSoilDiffusion
