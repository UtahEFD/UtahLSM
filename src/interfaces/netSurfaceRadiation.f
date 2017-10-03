      subroutine netSurfaceRadiation(surfaceTemperature,porosity,
     +     refTemp,surfaceMoisture,measRad,netRad,flag)
     
      use globals
      use SEBmodule
     
      integer*4 flag
      real*8 surfaceTemperature,refTemp,surfaceMoisture,netRad
      real*8,dimension(:) :: porosity,measRad
      end subroutine netSurfaceRadiation
