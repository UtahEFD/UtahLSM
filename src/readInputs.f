      subroutine readInputs(caseName)

      use globals
      use SEBmodule
      
      implicit none

      integer*4 i
      real*8 dtr,w,tkesgs,time
      character(len=32) :: caseName
      character(len=64) :: inputDir
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     READ in parameters from LESinputs !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      inputDir = 'inputs/' // caseName
                                          
      open(unit=1,file=trim(inputDir) // '/LSMinputs.txt',status='old')

      do i = 1,3
         read(1,*)
      enddo
      read(1,*) vonk
      read(1,*) pi
      do i=1,3
         read(1,*)
      enddo
      read(1,*) startUTC
      read(1,*) nsteps
      read(1,*) dtr
      read(1,*) z_m
      read(1,*) z_s
      read(1,*) z_i
      read(1,*) uScale
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
      do i = 1,3
         read(1,*)
      enddo

      ! READ soil type parameters
      read(1,*) soilLevels
	  read(1,*) z_o
      read(1,*) z_t
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
      read(1,*) albedo
      read(1,*) albedoMin
      do i = 1,3
         read(1,*) 
      enddo

      ! READ radiation parameters
      read(1,*) radiationFlag
      read(1,*) stepsPerRadVal
      read(1,*) SB_constant
      read(1,*) solarIrradiance
      read(1,*) lat
      read(1,*) long
      read(1,*) day
      read(1,*) emissivity

      close(1)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     READ in external soil and atmospheric data !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      allocate(zGnd(soilLevels),porosity(soilLevels),
     +     satPotential(soilLevels),satHydrCond(soilLevels),
     +     soilExponent(soilLevels),heatCapSoil(soilLevels),
     +     u(nsteps),v(nsteps))
      allocate(gndScalars(soilLevels,2),scalar(nsteps,scalarcount))

      ! soil
      write(*,*) 'reading in soil property data'
      Open (unit=1,file=trim(inputDir) // '/soilTypeParams.ini')
      read(1,*) porosity(1:soilLevels)
      read(1,*) satPotential(1:soilLevels)
      read(1,*) satHydrCond(1:soilLevels)
      read(1,*) soilExponent(1:soilLevels)
      read(1,*) heatCapSoil(1:soilLevels)
      close(1)
      
      write(*,*) 'reading in LSM parameter info'
      Open (unit=1,file=trim(inputDir) // '/soilLevels.ini')
      read(1,*) zGnd(1:soilLevels)
      close(1)

      write(*,*) 'reading in initial soil state'
      Open (unit=1,file=trim(inputDir) // '/soilTemperature.ini')
      read(1,*) gndScalars(1:soilLevels,1)
      close(1)
      Open (unit=1,file=trim(inputDir) // '/soilMoisture.ini')
      read(1,*) gndScalars(1:soilLevels,2)
      close(1)
      
      ! atmosphere    
      write(*,*) 'reading in external atmospheric data'      
      Open (unit=1,file=trim(inputDir) // '/timeseries_10.dat')
      do i=1,nsteps
         read(1,*) time,u(i),v(i),w,scalar(i,1),
     +        scalar(i,2),tkesgs
      enddo      

!!!!!!!!!!!!!!!!!!!!!!!!!
!     scalar parameters !
!!!!!!!!!!!!!!!!!!!!!!!!!

      ! scalar parameters, note these two could be assumed and 
      temperatureIndex=1
      moistureIndex=2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     nondimensionalization !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! time
      dt = dtr/(z_i/uScale)
      startUTC=startUTC*3600/(z_i/uScale)
      
      ! gravity
      g_hat=9.81d0*(z_i/(uScale**2))
      
      ! soil properties
      satPotential=satPotential/z_i
      satHydrCond=satHydrCond/uScale
      heatCapSoil=heatCapSoil/(densityAir*Cp_air)
      zGnd=zGnd/z_i
      
      ! atm length scales
      z_o=z_o/z_i
      z_t=z_t/z_i 
      z_m=z_m/z_i
      z_s=z_s/z_s
      
      ! atm velocity and scalars
      u=u/uScale
      v=v/uScale
      do i=1,scalarcount
         scalar(:,i)=scalar(:,i)/scalarScales(i)
      enddo
      
      ! water properties
      densityWater = densityWater/densityAir
      latentHeatWater = 
     +     latentHeatWater/(Cp_air*scalarScales(temperatureIndex))
      heatCapWater = heatCapWater/(densityAir*Cp_air)
      waterGasConst = 
     +     waterGasConst*scalarScales(temperatureIndex)/uScale**2
      
      ! radiation model
      SB_constant=SB_constant*scalarScales(temperatureIndex)**3/
     +     (Cp_air*densityAir*uScale)
      solarIrradiance = 
     +     solarIrradiance/(scalarScales(temperatureIndex)*uScale)
      lat = lat*pi/180.d0
      long = long*pi/180.d0   
      
      return
      end
