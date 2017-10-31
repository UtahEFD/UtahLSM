!     LSMmodules.f contains all the required modules for utahLSM.f
!     Each module declares the variables required for a specific set of 
!     All subroutines use globals
!     The modules in this file include: 
!              globals
!              SEBmodule
!
      MODULE globals
      implicit none
      real*8 dt
      integer*4 t
      integer*4 nsteps,scalarCount
      real*8 Pi,grav,vonk,z_i, z_m, z_s
      real*8 UTC,startUTC,UTC_hrs
      integer*4 temperatureIndex,moistureIndex
      real*8,allocatable,dimension(:)::surfaceFluxes,scalarScales,
     +     zGnd,porosity,satPotential,satHydrCond,soilExponent,
     +     heatCapSoil,u,v,measRad
      real*8,allocatable,dimension(:,:) :: gndScalars,scalar
      
      END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      MODULE SEBmodule

      integer*4 soilLevels,maxFluxIterations,maxTempIterations,
     +     endConstSEB,updateFreqSEB,integrateSoilDiffFreq,
     +     radiationFlag,stepsPerRadVal,day,albedoFlag,
     +     tempConvFlag,bisectFlag,sepFlag,computeLH
      real*8 z_o,z_t,pressureScale,densityAir,Cp_air,densityWater,
     +     latentHeatWater,heatCapWater,waterGasConst,moistureCriteria,
     +     temperatureCriteria,tempFluxCriteria,convFactor,SB_constant,
     +     solarIrradiance,lat,long,emissivity,albedo,albedoMin

      END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!