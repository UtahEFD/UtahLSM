!     LSMmodules.f contains all the required modules for LSM_ShingletonV1.f
!     Each module declares the variables required for a specific set of 
!     All subroutines use globals
!     The modules in this file include: 
!              globals
!              SEBmodule
!
      MODULE globals
!     ALL program units will 'use globals'
      implicit none

      real*8 dt
      integer*4 t
      integer*4 nsteps,scalarCount
      real*8 Pi,g_hat,u_star,vonk,z_i
      real*8 UTC,startUTC,UTC_hrs
      integer*4 temperatureIndex,moistureIndex
      real*8,allocatable,dimension(:)::surfaceFluxes,scalarScales
      
      END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      MODULE SEBmodule

      integer*4 soilLevels,maxFluxIterations,maxTempIterations,
     +     endConstSEB,updateFreqSEB,integrateSoilDiffFreq,
     +     radiationFlag,stepsPerRadVal,day,albedoFlag
      real*8 zt,pressureScale,densityAir,Cp_air,densityWater,
     +     latentHeatWater,heatCapWater,waterGasConst,moistureCriteria,
     +     temperatureCriteria,tempFluxCriteria,convFactor,SB_constant,
     +     solarIrradiance,lat,long,emissivity

      END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
