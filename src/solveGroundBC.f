      subroutine solveGroundBC(Uref,scalarRef, 
     +     gndScalars,ustar,scalarFlux,soilHeatFlux,porosity,
     +     satPotential,satHydrCond,soilExponent,heatCapSoil,
     +     measRad,netRad,Psi,Psi0,fi,fiH,zGnd)
      
      use globals
      use SEBmodule
      implicit none
      
      interface
         include './interfaces/getStabilityCorrections.f'
         include './interfaces/getSurfaceMixingRatio.f'
         include './interfaces/integrateSoilDiffusion.f'
         include './interfaces/getSoilThermalTransfer.f'
         include './interfaces/getWaterConductivity.f'
         include './interfaces/netSurfaceRadiation.f'
      end interface

      real*8 Uref,ustar,soilHeatFlux,
     +     netRad,Psi,Psi0,fi,fiH,PsiH,PsiH0
      real*8,dimension(:) :: scalarRef,scalarFlux,measRad,
     +     zGnd,porosity,satPotential,satHydrCond,soilExponent,
     +     heatCapSoil
      real*8,dimension(:,:) :: gndScalars

      integer*4 bisectFlag,iterateFlux,iterateTemp,
     +     i,k,tempConvFlag ,sepFlag,computeLH,
     +     skipSEBFlag,skipIntegrationFlag
      real*8, dimension(2)::moistPotential,
     +     conductivity,diffusConduct,hydrConduct,lastSurfScalars,
     +     lastScalarFlux
      real*8 denom,denomH,x y,z,deta_dz,h,partialPressure,qs,q_gnd,
     +     specHum,specHum_gnd,shortIn,shortNet,longNet,longIn,longOut,
     +     soilMoistureFlux,mFluxConvergence,TsConvergence,
     +     tFluxConvergence,measResidual,measPress,lastTemperature,
     +     lastSoilMoistureFlux,SEB,SEBlast,dSEB_dT,lastSoilFlux,
     +     lastTempFlux,lowT,highT,dTs_Old,dTs,dFluxM_dT, dFluxT_dT
      real*8, dimension(size(gndScalars,1)) :: tempT, tempM
      
!!! These flags control some of the iteration and methods (example computeLH turns off
!!! the latent heat calc

      tempConvFlag = 1
      bisectFlag   = 0
      sepFlag      = 1
      computeLH    = 1

      lastSurfScalars(:) = gndScalars(1,:)
      lastTemperature    = gndScalars(1,temperatureIndex)

      if ( UTC.eq.(startUTC + dt).or.t==1 )then ! first timeStep
         lowT    = gndScalars(1,temperatureIndex) - 5/scalarScales(1) ! set bounds on temperature 
         highT   = gndScalars(1,temperatureIndex) + 5/scalarScales(1) ! for convergence
      else
!     if solution converged 1st time step, it shouldn't change much there after
         lowT    = gndScalars(1,temperatureIndex) - 0.25/scalarScales(1) ! set bounds on temperature 
         highT   = gndScalars(1,temperatureIndex) + 0.25/scalarScales(1) ! for convergence
      endif
      dTs_Old = highT - lowT
      dTs     = dTs_Old

      skipSEBFlag = 0
      if( t > endConstSEB .and.
     >     mod(t-endConstSEB,updateFreqSEB) /= 0 )then
         skipSEBFlag = 1
      endif

      if(skipSEBFlag.eq.1)then
         write(*,*)'skipping SEB, t=',t
      endif
      
      if ( UTC == (startUTC + dt).or.t==1 .or. skipSEBFlag == 1 )then

         Psi=0.d0
         Psi0=0.d0
         PsiH=0.d0
         PsiH0=0.d0

         do i=1,4

!           momentum
            denom = dlog( z_m / zo ) + Psi - Psi0
            ustar = Uref*vonk/denom
            
!           scalar flux
            denomH = dlog( z_s / zt ) + PsiH - PsiH0 
                   
            scalarFlux(temperatureIndex) = 
     >           ( gndScalars( 1,temperatureIndex )
     >           - scalarRef(temperatureIndex) ) * ustar*vonk/denomH
            
            call getSurfaceMixingRatio(gndScalars,q_gnd,measPress,
     >           porosity,satPotential,soilExponent)    
            if(computeLH == 1) then
               scalarFlux(moistureIndex) = ( q_gnd 
     >              - scalarRef(moistureIndex) )*ustar*vonk/denomH
               
            endif

! compute Psi and fi values for momentum and scalars from computed flux

            call getStabilityCorrections(ustar,
     >           scalarRef(temperatureIndex),
     >           scalarRef(moistureIndex),scalarFlux,
     >           Psi,Psi0,fi,PsiH,PsiH0,fiH)
            
         enddo

      endif   
      
! solve ground budget iteratively for temperature and moisture
!     with solution of PsiH and PsiH0, compute flux of other scalars (if any). 
      
      if(skipSEBFlag == 0)then
         do iterateFlux = 1,maxFluxIterations
            if( computeLH==1 ) then
!     assume surface fluxes are constant and converge to surface temperature
!     using Newton-Raphson, where SEB = f(Ts)         
               do iterateTemp = 1, maxTempIterations

!     SURFACE ENERGY BUDGET
!     solve surface energy budget 
                  
                  call netSurfaceRadiation(gndScalars(1,
     >                 temperatureIndex),porosity,
     >                 scalarRef(temperatureIndex),
     >                 gndScalars(1,moistureIndex),measRad,
     >                 netRad,iterateFlux*iterateTemp)

                  call getSoilThermalTransfer(
     >                 gndScalars(1:2,moistureIndex),
     >                 conductivity,0,porosity(1:2),satPotential(1:2),
     >                 soilExponent(1:2),heatCapSoil(1:2))

                  soilHeatFlux = ( gndScalars(1,temperatureIndex) - 
     >                 gndScalars(2,temperatureIndex) )*
     >                 (sum(conductivity)/2.d0)/ (zGnd(2) - zGnd(1))
                  
                  SEB = soilHeatFlux + scalarFlux(temperatureIndex)
     >                 + scalarFlux(moistureIndex)*latentHeatWater
     >                 - netRad !- shortNet - longNet
!     netRad = shortNet + longNet
                 
!     compute dSEB_dT              
                  dFluxT_dT = ustar*vonk/denomH 
                  
                  dSEB_dT = 4.d0*emissivity*SB_Constant*
     >                 gndScalars(1,temperatureIndex)**3
     >                 + (sum(conductivity)/2)/zGnd(2) + dFluxT_dT
                  
                  dTs = SEB / dSEB_dT
                  
                  gndScalars(1,temperatureIndex) = 
     >                 gndScalars(1,temperatureIndex) - dTs

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
     >           - scalarRef(temperatureIndex) ) * ustar*vonk/denomH
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
            
            do i=1,4

!              momentum 
               denom = dlog( z_m / zo ) + Psi - Psi0
               ustar = Uref*vonk/denom
               
!              scalar flux
               denomH = dlog( z_s / zt ) + PsiH - PsiH0               
               
               scalarFlux(temperatureIndex) = ( 
     >              gndScalars(1,temperatureIndex)
     >              - scalarRef(temperatureIndex) ) * ustar*vonk/denomH
         
               call getSurfaceMixingRatio(gndScalars,q_gnd,measPress,
     >              porosity,satPotential,soilExponent)      
               if(computeLH == 1 )then
                  scalarFlux(moistureIndex) = ( q_gnd 
     >                 - scalarRef(moistureIndex) )*ustar*vonk/denomH
               endif
         
!     compute Psi and fi values for momentum and scalars from computed flux

               call getStabilityCorrections(ustar,
     >              scalarRef(temperatureIndex),scalarRef(moistureIndex)
     >              ,scalarFlux,Psi,Psi0,fi,PsiH,PsiH0,fiH) 
           
            enddo
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

                  call getSurfaceMixingRatio(gndScalars,q_gnd,measPress,
     >                 porosity,satPotential,soilExponent)      
                  if(computeLH == 1 )then
                     scalarFlux(moistureIndex) = ( q_gnd 
     >                    - scalarRef(moistureIndex) )*ustar*vonk/denomH
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
            call integrateSoilDiffusion(gndScalars,
     >           lastSurfScalars,1,zGnd,porosity,satPotential,
     >           soilExponent,heatCapSoil,satHydrCond)
         endif
         
         if (moistureIndex /= 0)then
            call integrateSoilDiffusion(gndScalars,
     >           lastSurfScalars,2,zGnd,porosity,satPotential,
     >           soilExponent,heatCapSoil,satHydrCond)
         endif
      endif
      
      end subroutine solveGroundBC
