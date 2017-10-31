      subroutine solveGroundBC(Uref,scalarRef,ustar,scalarFlux,
     +     soilHeatFlux,netRad)
      
      use globals
      use SEBmodule
      implicit none
      
      interface
         include './interfaces/getFluxesMOST.f'
         include './interfaces/getSurfaceMixingRatio.f'
         include './interfaces/integrateSoilDiffusion.f'
         include './interfaces/getSoilThermalTransfer.f'
         include './interfaces/getWaterConductivity.f'
         include './interfaces/netSurfaceRadiation.f'
      end interface

      real*8 Uref,ustar,soilHeatFlux,netRad,psi,psi0,fi,fiH,psiH,psiH0
      real*8,dimension(:) :: scalarRef,scalarFlux
      
      integer*4 iterateFlux,iterateTemp,
     +     i,k,skipSEBFlag,skipIntegrationFlag
      real*8, dimension(2)::moistPotential,
     +     conductivity,diffusConduct,hydrConduct,lastSurfScalars,
     +     lastScalarFlux
      real*8 x y,z,deta_dz,h,partialPressure,qs,q_gnd,
     +     specHum,specHum_gnd,shortIn,shortNet,longNet,longIn,longOut,
     +     soilMoistureFlux,mFluxConvergence,TsConvergence,
     +     tFluxConvergence,measResidual,measPress,lastTemperature,
     +     lastSoilMoistureFlux,SEB,SEBlast,dSEB_dT,lastSoilFlux,
     +     lastTempFlux,lowT,highT,dTs_Old,dTs,dFluxM_dT, dFluxT_dT
      real*8, dimension(size(gndScalars,1)) :: tempT, tempM
      
      ! set convergence bounds on temperature 
      ! if solution converged on 1st time step, it shouldn't change 
      ! much thereafter
      if ( UTC.eq.(startUTC + dt).or.t==1 )then
         lowT    = gndScalars(1,temperatureIndex) - 5
         highT   = gndScalars(1,temperatureIndex) + 5
      else
         lowT    = gndScalars(1,temperatureIndex) - 0.25 
         highT   = gndScalars(1,temperatureIndex) + 0.25
      endif
      
      dTs_Old = highT - lowT
      dTs     = dTs_Old
      
      ! store previous scalars for convergence tests
      lastSurfScalars(:) = gndScalars(1,:)
      lastTemperature    = gndScalars(1,temperatureIndex)
      
      ! determine whether to compute SEB 
      skipSEBFlag = 0
      
      if( t > endConstSEB .and.
     >     mod(t-endConstSEB,updateFreqSEB) /= 0 )then
         skipSEBFlag = 1
         write(*,*)'skipping SEB, t=',t
      endif
      
      ! if first time step or in-between SEB calls, just use MOST
      if ( UTC == (startUTC + dt).or.t==1 .or. skipSEBFlag == 1 )then

         psi   = 0.d0
         psi0  = 0.d0
         psiH  = 0.d0
         psiH0 = 0.d0
         
         call getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
     >        psi,psi0,fi,psiH,psiH0,fiH)

      endif
      
      ! solve SEB iteratively for temperature and moisture
      ! with solution of psiH and psiH0, compute flux of other scalars  
      if(skipSEBFlag == 0) then
         
         ! begin: outer flux iteration
         do iterateFlux = 1,maxFluxIterations
            
            ! begin: inner temperature iteration
            ! assume surface fluxes are constant and converge to Ts
            ! using Newton-Raphson, where SEB = f(Ts)         
            do iterateTemp = 1, maxTempIterations

               ! compute net surface radiation                   
               call netSurfaceRadiation(scalarRef(temperatureIndex),
     >              gndScalars(1,moistureIndex),netRad)
               
               ! compute soil conductivity K
               call getSoilThermalTransfer(0,
     >              gndScalars(1:2,moistureIndex),
     >              porosity(1:2),satPotential(1:2),
     >              soilExponent(1:2),heatCapSoil(1:2),conductivity)
                                 
               ! compute soil heat flux
               soilHeatFlux = -( gndScalars(2,temperatureIndex) - 
     >              gndScalars(1,temperatureIndex) )*
     >              (sum(conductivity)/2.d0)/ (zGnd(2) - zGnd(1))

               ! compute surface energy balance
               SEB = soilHeatFlux
     >              + scalarFlux(temperatureIndex)*densityAir*Cp_air
     >              + scalarFlux(moistureIndex)
     >              * densityAir*latentHeatWater - netRad

               ! compute derivative of heat flux              
               dFluxT_dT = (ustar*vonk / ( dlog( z_s / z_t )
     >                   - psiH + psiH0 ))*densityAir*Cp_air
                 
               ! compute time derivative of SEB
               dSEB_dT = 4.d0*emissivity*SB_Constant*
     >              gndScalars(1,temperatureIndex)**3
     >              + (sum(conductivity)/2)/zGnd(2) + dFluxT_dT
               
               ! compute time change in surface temperature
               dTs = SEB / dSEB_dT
               
               ! update surface temperature
               gndScalars(1,temperatureIndex) = 
     >              gndScalars(1,temperatureIndex) - dTs
               
               ! check for convergence
               TsConvergence = 
     >              abs(dTs)/gndScalars(1,temperatureIndex)
               if ( TsConvergence < temperatureCriteria )then
                  print*,"New Temp found:",iterateTemp,
     >                 gndScalars(1,temperatureIndex)
                  exit
               endif
            
            ! end: inner temperature iteration 
            enddo                  
            
            ! temperature converged, now solve heat flux          
            ! save previous flux for convergence test
            lastTempFlux = scalarFlux(temperatureIndex)

            ! re-compute heat flux using new temperature
!            call getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
!     >           psi,psi0,fi,psiH,psiH0,fiH)
            scalarFlux(temperatureIndex) = ( 
     >           gndScalars(1,temperatureIndex )
     >           - scalarRef(temperatureIndex) ) 
     >           * ustar*vonk/(dlog( z_s / z_t ) - psiH + psiH0 )
            
            ! check for convergence
            tFluxConvergence = 
     >           abs( scalarFlux(temperatureIndex)-lastTempFlux)
            if (((tFluxConvergence < tempFluxCriteria)).or.
     >           (iterateFlux >= maxFluxIterations)) then
               print*, "New Flux: ",scalarFlux(temperatureIndex)
               exit
            endif
            
            ! if no convergence, bisect temperature
            gndScalars(1,temperatureIndex) = (lastTemperature+
     >           gndScalars(1,temperatureIndex))/2.d0
            
            ! save bisected temperature for next loop
            lastTemperature = gndScalars(1,temperatureIndex)
            
            ! re-compute new fluxes and stability functions
            call getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
     >           psi,psi0,fi,psiH,psiH0,fiH)
         
         ! end: outer flux iteration  
         enddo
      
      ! end: if computing SEB
      endif
!         
!         ! heat flux converged, hold constant, compute latent heat flux
!         ! compute water conductivity
!         call getWaterConductivity(gndScalars(1:2,moistureIndex),
!     >        porosity(1:2),satPotential(1:2),satHydrCond(1:2),
!     >        soilExponent(1:2),diffusConduct,hydrConduct)
!         
!         ! compute soil moisture flux
!         soilMoistureFlux = densityWater*(sum(diffusConduct)/2.d0)
!     >        *(gndScalars(1,moistureIndex)-gndScalars(2,moistureIndex))
!     >        /(zGnd(2)-zGnd(1)) + densityWater*sum(hydrConduct)/2.d0
!
!         ! convergence loop
!         do i = 1,maxTempIterations
!            
!            ! get surface mixing ratio
!            call getSurfaceMixingRatio(q_gnd)      
!            
!            ! compute latent heat flux
!            scalarFlux(moistureIndex) = ( q_gnd 
!     >           - scalarRef(moistureIndex) )* ustar
!     >           * vonk/( dlog( z_s / z_t ) - psiH + psiH0 )
!            
!            ! compute water conductivity
!            call getWaterConductivity(
!     >           gndScalars(1:2,moistureIndex),
!     >           diffusConduct,hydrConduct,porosity(1:2),
!     >           satPotential(1:2),satHydrCond(1:2),
!     >           soilExponent(1:2))
!            
!            ! save soil moisture flux
!            lastSoilMoistureFlux = soilMoistureFlux
!            
!            ! compute weighted soil moisture flux
!            soilMoistureFlux = convFactor*lastSoilMoistureFlux - 
!     >           (1.d0-convFactor)*scalarFlux(moistureIndex)
!            
!            ! compute for new surface moisture potential
!            moistPotential(2) = satPotential(2)*(porosity(2)
!     >           /gndScalars(2,moistureIndex))**soilExponent(2)
!            
!            moistPotential(1) = moistPotential(2) + 
!     >           (zGnd(2)-zGnd(1))
!     >           *((soilMoistureFlux/
!     >           (densityWater*sum(hydrConduct)/2.d0))-1.d0)
!            
!            ! compute new surface moisture content
!            gndScalars(1,moistureIndex) = porosity(1)*
!     >           (moistPotential(1)/satPotential(1))
!     >           **(-1.d0/soilExponent(1))
!            
!            ! check for moisture convergence
!            mFluxConvergence = abs(
!     >           (scalarFlux(moistureIndex)
!     >           - (soilMoistureFlux))
!     >           /(scalarFlux(moistureIndex)))  
!            if( mFluxConvergence < moistureCriteria )then
!               print*, "Fuck yeah"
!               exit
!            endif
!         enddo
!      endif
!
!      skipIntegrationFlag = 0
!      if( t > endConstSEB .and.
!     >     mod(t-endConstSEB,integrateSoilDiffFreq) /= 0 )then
!         skipIntegrationFlag = 1
!      endif
!      
!      ! use surface fluxes and soil profiles to integrate 
!      ! heat diffusion equation and soil moisture equation in time
!      if( skipIntegrationFlag == 0 )then
!         if (temperatureIndex /= 0)then
!            call integrateSoilDiffusion(lastSurfScalars,1)
!         endif
!         
!         if (moistureIndex /= 0)then
!            call integrateSoilDiffusion(lastSurfScalars,2)
!         endif
!      endif
      
      end subroutine solveGroundBC