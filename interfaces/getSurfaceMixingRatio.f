      subroutine getSurfaceMixingRatio(gndScalars,q_gnd,measPress,
     +     porosity,satPotential,soilExponent)
      real*8 q_gnd,measPress
      real*8,dimension(:) :: porosity,satPotential,soilExponent
      real*8,dimension(:,:) :: gndScalars
      end subroutine getSurfaceMixingRatio
