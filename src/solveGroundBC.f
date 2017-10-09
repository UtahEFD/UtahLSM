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

      real*8 Uref,ustar,soilHeatFlux,netRad,Psi,Psi0,fi,fiH,PsiH,PsiH0
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
      ! if solution converged 1st time step, it shouldn't change much there after
      if ( UTC.eq.(startUTC + dt).or.t==1 )then
         lowT    = gndScalars(1,temperatureIndex) - 5/scalarScales(1)
         highT   = gndScalars(1,temperatureIndex) + 5/scalarScales(1)
      else
         lowT    = gndScalars(1,temperatureIndex) - 0.25/scalarScales(1) 
         highT   = gndScalars(1,temperatureIndex) + 0.25/scalarScales(1)
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

         Psi=0.d0
         Psi0=0.d0
         PsiH=0.d0
         PsiH0=0.d0
         
         call getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
     >        Psi,Psi0,fi,PsiH,PsiH0,fiH)

      endif   
      
      ! solve SEB iteratively for temperature and moisture
      ! with solution of PsiH and PsiH0, compute flux of other scalars (if any).  
      if(skipSEBFlag == 0)then
         do iterateFlux = 1,maxFluxIterations
            if( computeLH==1 ) then
               ! assume surface fluxes are constant and converge to surface temperature
               ! using Newton-Raphson, where SEB = f(Ts)         
               do iterateTemp = 1, maxTempIterations

                  ! solve surface energy budget                   
                  call netSurfaceRadiation(scalarRef(temperatureIndex),
     >                 netRad)

                  call getSoilThermalTransfer(
     >                 gndScalars(1:2,moistureIndex),
     >                 conductivity,0,porosity(1:2),satPotential(1:2),
     >                 soilExponent(1:2),heatCapSoil(1:2))

                  soilHeatFlux = ( gndScalars(1,temperatureIndex) - 
     >                 gndScalars(2,temperatureIndex) )*
     >                 (sum(conductivity)/2.d0)/ (zGnd(2) - zGnd(1))
                  
                  ! netRad = shortNet + longNet
                  SEB = soilHeatFlux + scalarFlux(temperatureIndex)
     >                 + scalarFlux(moistureIndex)*latentHeatWater
     >                 - netRad
                 
                  ! compute dSEB_dT              
                  dFluxT_dT = ustar*vonk / ( dlog( z_s / z_t )
     >                      + PsiH - PsiH0 )
                  
                  dSEB_dT = 4.d0*emissivity*SB_Constant*
     >                 gndScalars(1,temperatureIndex)**3
     >                 + (sum(conductivity)/2)/zGnd(2) + dFluxT_dT
                  
                  dTs = SEB / dSEB_dT
                  
                  gndScalars(1,temperatureIndex) = 
     >                 gndScalars(1,temperatureIndex) - dTs
                  
                  if (t<2) then 
                     print*,SEB
                  else    
                     call abort
                  endif 

                  TsConvergence = 
     >                 abs(dTs)/gndScalars(1,temperatureIndex)
                  
                  if ( TsConvergence < temperatureCriteria )then
                     exit
                  endif
            
                  SEBlast = SEB        
                  lastScalarFlux = scalarFlux
                  lastSoilFlux   = soilHeatFlux
               enddo            ! iterateTemp                  
             
            else
            
               call getSoilThermalTransfer(gndScalars(1:2,moistureIndex)
     >              ,conductivity,0,porosity(1:2),satPotential(1:2),
     >              soilExponent(1:2),heatCapSoil(1:2))
               
               soilHeatFlux = ( gndScalars(1,temperatureIndex) - 
     >              gndScalars(2,temperatureIndex) )*(sum(conductivity)/
     >              2.d0)/(zGnd(2) - zGnd(1))
            endif
            
            lastTempFlux = scalarFlux(temperatureIndex)
            
            scalarFlux(temperatureIndex) = ( 
     >           gndScalars( 1,temperatureIndex )
     >           - scalarRef(temperatureIndex) ) 
     >           * ustar*vonk/(dlog( z_s / z_t ) + PsiH - PsiH0 )
            tFluxConvergence = 
     >           abs( scalarFlux(temperatureIndex)-lastTempFlux)
            
            if( sepFlag == 1 )then
               mFluxConvergence = 0
            endif

            if( ((mFluxConvergence < moistureCriteria) .and.
     >           (tFluxConvergence < tempFluxCriteria)).or.
     >           (iterateFlux >= maxFluxIterations) )then
               exit
            endif
            
            gndScalars(1,temperatureIndex) = (lastTemperature+
     >           gndScalars(1,temperatureIndex))/2.d0
            
            lastTemperature = gndScalars(1,temperatureIndex)
            
            call getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
     >           Psi,Psi0,fi,PsiH,PsiH0,fiH)
           
         enddo                  ! iterate=1,maxFluxIterations   
 
!    compute first soil moisture flux
         call getWaterConductivity(gndScalars(1:2,moistureIndex),
     >        diffusConduct,hydrConduct,porosity(1:2),satPotential(1:2),
     >        satHydrCond(1:2),soilExponent(1:2))
         
         soilMoistureFlux = densityWater*(sum(diffusConduct)/2.d0)
     >        *(gndScalars(1,moistureIndex)-gndScalars(2,moistureIndex))
     >        /(zGnd(2)-zGnd(1)) + densityWater*sum(hydrConduct)/2.d0
         
         if( computeLH == 1 )then
            if( sepFlag == 1 )then
               do i = 1,200     !maxTempIterations

                  call getSurfaceMixingRatio(q_gnd)      
                  if(computeLH == 1 )then
                     scalarFlux(moistureIndex) = ( q_gnd 
     >                    - scalarRef(moistureIndex) )* ustar
     >                    * vonk/( dlog( z_s / z_t ) + PsiH - PsiH0 )
                  endif

                  call getWaterConductivity(
     >                 gndScalars(1:2,moistureIndex),
     >                 diffusConduct,hydrConduct,porosity(1:2),
     >                 satPotential(1:2),satHydrCond(1:2),
     >                 soilExponent(1:2))
            
                  lastSoilMoistureFlux = soilMoistureFlux
                  soilMoistureFlux = convFactor*lastSoilMoistureFlux - 
     >                 (1.d0-convFactor)*scalarFlux(moistureIndex)
            
                  moistPotential(2) = satPotential(2)*(porosity(2)
     >                 /gndScalars(2,moistureIndex))**soilExponent(2)
            
                  moistPotential(1) = moistPotential(2) + 
     >                 (zGnd(2)-zGnd(1))
     >                 *((soilMoistureFlux/
     >                 (densityWater*sum(hydrConduct)/2.d0))-1.d0)
                  
                  gndScalars(1,moistureIndex) = porosity(1)*
     >                 (moistPotential(1)/satPotential(1))
     >                 **(-1.d0/soilExponent(1))

                  mFluxConvergence = abs(
     >                 (scalarFlux(moistureIndex)
     >                 - (-soilMoistureFlux))
     >                 /(scalarFlux(moistureIndex)))
                  
                  if( mFluxConvergence < moistureCriteria )then
                     exit
                  endif
               enddo
            endif
         endif
      endif                     ! if( skipSEBFlag == 0 )

      skipIntegrationFlag = 0
      if( t > endConstSEB .and.
     >     mod(t-endConstSEB,integrateSoilDiffFreq) /= 0 )then
         skipIntegrationFlag = 1
      endif
      
!     use surface fluxes and soil profiles to integrate 
!     heat diffusion equation and soil moisture equation in time
      if( skipIntegrationFlag == 0 )then
         if (temperatureIndex /= 0)then
            call integrateSoilDiffusion(lastSurfScalars,1)
         endif
         
         if (moistureIndex /= 0)then
            call integrateSoilDiffusion(lastSurfScalars,2)
         endif
      endif
      
      end subroutine solveGroundBC
