      subroutine netSurfaceRadiation(surfaceTemperature,albedo,
     +     minAlbedo,porosity,refTemp,surfaceMoisture,
     +     measRad,netRad,flag)
      integer*4 flag
      real*8 surfaceTemperature,refTemp,surfaceMoisture,netRad,
     +     albedo,minAlbedo
      real*8,dimension(:) :: porosity,measRad
      end subroutine netSurfaceRadiation
