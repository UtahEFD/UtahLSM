      subroutine netSurfaceRadiation(refTemp,sfcMois,netRad)
      ! module: computes the net surface radiation
      !         = incoming lw - outgoing lw + incoming sw - outgoing sw
      !         
      !         incoming sw - entirely due to solar radiation
      !         outgoing sw - due to reflection of solar radiation
      !                       dependent on surface albedo
      !         incoming lw - due to radiation emitted by atm/env
      !                       currently neglected
      !         outgoing lw - due to emittence of the surface
      !                       computed using Stefan-Boltzmann
      !
      !         supports the use of measured net radiation
      !
      ! inputs: refTemp - measured near-surface temerature
      !         sfcMois - measured surface moisture
      !         measRad - external net radiation (radiationFlag=1)
      !         
      ! output: netRad  - net radiation
      
      use globals
      use SEBmodule
      implicit none

      interface
         include './interfaces/solarRadiation.f'
      end interface

      integer*4 flag
      real*8 refTemp,sfcMois,netRad
      
      integer*4 i,j
      real*8  shortNet, longNet
     
      ! if using measured radiation, interpolate and return
      if( radiationFlag == 1 ) then

         i = mod(t-1,stepsPerRadVal)
         j = int( (t - i)/stepsPerRadVal ) + 1

         netRad = measRad(j) + i*( measRad(j+1) - measRad(j) )
     >        / stepsPerRadVal
         return
      endif
      
      ! estimate for net longwave
      longNet = -0.04d0
      
      ! compute net shortwave radiation
      call solarRadiation(shortNet,sfcMois)
      
      ! compute net radiation
      netRad = shortNet + longNet
            
      end subroutine netSurfaceRadiation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine solarRadiation(shortNet,sfcMois)
      ! solarRadiation approximates the incoming solar radiation 
      ! based on the latitude, longitude, day of year and time of day
      ! neglects attenuation due to cloud cover (ie assumes clear sky)
      ! neglects attenuation due to moist air or air pollution
      ! assumes slope of ground surface is zero and 
      ! earths rotation around the sun is circular
      !
      ! inputs:
      !     lat  - latitude of location to compute radiation (degrees)
      !     long - longitude of location to compute radiation (degrees)
      !     day  - julian day of year (1-365)
      !     UTC  - Coordinated Universal Time (UTC) of day (hours)
      ! outputs:
      !     radiation - solar short wave radiation flux at the
      !                 earths surface given the above inputs (negative is upward, positive is downward).
      ! 
      ! model reference: Stull, Roland B., An Introduction to Boundary Layer Meteorology.
      !                         pg. 257
      !
      use globals
      use SEBmodule
      implicit none

      real*8 shortNet,sfcMois

      real*8 Az, As, A, E, Z, declination, sinElevation,
     +     transmissivity

      declination = 23.45d0*(pi/180.d0)*cos(2.d0*pi*(day-173)/365.25d0)

      sinElevation = sin(lat)*sin(declination) - 
     >     cos(lat)*cos(declination)*
     >     cos( (2*pi*UTC/(24.d0*3600.d0)) - long )

      if(sinElevation > 0)then
         transmissivity = (0.6d0 + 0.2d0*sinElevation)
          ! compute albedo based on material, moisture content and sinElevation
          ! ADD computation for albedo based on moisture and angle
         if( albedoFlag == 1 )then
            E  = asin(sinElevation)            
            Z  = (pi/2.d0) - E
            Az = 0.01d0*(exp(0.003286d0*Z**1.5d0) - 1.d0)
            if( sfcMois/porosity(1) <= 0.5)then
               As = albedo*(1 - sfcMois/porosity(1))
            else
               As = albedoMin
            endif
            A  = As + Az
            ! for gables3 (from PILPS paper, online say albedo = 0.23 ??)
            ! A = albedo - albedoMin*sinElevation 
         else
            A = albedo
         endif

         shortNet = (1.0-A)*solarIrradiance*transmissivity*sinElevation
      else
         shortNet = 0
      endif

      end subroutine solarRadiation