      subroutine readInputs()

      use globals
      use SEBmodule
      
      implicit none

      integer*4 i
      real*8 dtr

!     READ in parameters 
                                                       
      open(unit=1,file='inputs/LSMinputs.txt',status='old')

      do i = 1,6
         read(1,*)
      enddo
      read(1,*) vonk
      read(1,*) pi
      do i=1,3
         read(1,*)
      enddo
      read(1,*) nsteps
      read(1,*) dtr
      read(1,*) z_i
      read(1,*) u_star
      do i = 1,3
         read(1,*)
      enddo

! READ scalar parameters               
      read(1,*) scalarCount
      if(scalarCount.gt.0)then
      ALLOCATE( scalarScales(scalarCount))
      do i=1,scalarCount
         read(1,*) scalarScales(i)
      enddo
      endif
!     divide by scalarScales for temperature                            
      do i = 1,3
         read(1,*)
      enddo

!     READ soil type parameters

      read(1,*) soilLevels

      read(1,*) zt
      read(1,*) pressureScale
      read(1,*) densityAir
      read(1,*) Cp_air
      read(1,*) densityWater
      read(1,*) latentHeatWater
      read(1,*) heatCapWater
      read(1,*) waterGasConst
      read(1,*) moistureCriteria
      read(1,*) temperatureCriteria
      read(1,*) tempFluxCriteria
      read(1,*) maxFluxIterations
      read(1,*) maxTempIterations
      read(1,*) convFactor
      read(1,*) endConstSEB
      read(1,*) updateFreqSEB
      read(1,*) integrateSoilDiffFreq
      read(1,*) albedoFlag
      do i = 1,3
         read(1,*) 
      enddo

!     READ radiation parameters

      read(1,*) radiationFlag
      read(1,*) stepsPerRadVal
      read(1,*) SB_constant
      read(1,*) solarIrradiance
      read(1,*) lat
      read(1,*) long
      read(1,*) day
      read(1,*) emissivity

      close(1)

!     calculated grid parameters
      dt = dtr/(z_i/u_star)
      startUTC = 0.d0*3600/(z_i/u_star)

!     scalar parameters, note these two could be assumed and 
   
      temperatureIndex=1
      moistureIndex=2

!     nondimensionalizations

      g_hat=9.81d0*(z_i/(u_star**2))

      if(soilLevels.gt.0)then
      densityWater = densityWater/densityAir
      latentHeatWater = 
     +     latentHeatWater/Cp_air*scalarScales(temperatureIndex)
      heatCapWater = heatCapWater/(densityAir*Cp_air)
      waterGasConst = 
     +     waterGasConst*scalarScales(temperatureIndex)/u_star**2
      SB_constant=SB_constant*scalarScales(temperatureIndex)**3/
     +     (Cp_air*densityAir*u_star)
      solarIrradiance = 
     +     solarIrradiance/(scalarScales(temperatureIndex)*u_star)
      lat = lat*pi/180.d0
      long = long*pi/180.d0   
      endif

      return
      end
