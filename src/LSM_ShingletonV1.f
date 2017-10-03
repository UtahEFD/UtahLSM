c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     c    This program is developed to run the LSM          c
c     c    model developed by Shingleton (2010)              c
c     c                                                      c
c     c                 October, 2012                        c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Program LSM_ShingletonV1

      use globals
      use SEBmodule
      implicit none
      
      interface
         include './interfaces/solveGroundBC.f'
      end interface

!!!!!!Variable Declarations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c short integer and character variables

      integer*4 i,j
      
c double precision scalars
      real*8 M,w,tkesgs,time,Psi,
     +     Psi0,fi,fi_H,ustar,soilHeatFlux,netRad

c double precision 1D vectors
      real*8,allocatable,dimension(:) :: zGnd,porosity,satPotential,
     +     satHydrCond,soilExponent,heatCapSoil,u,v,
     +     measRad,scalarFlux

c double precision 2D arrays
      real*8,allocatable,dimension(:,:) :: gndScalars,scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'reading inputs from LSMinputs file'
      call readInputs()

!!!!!!!!!!!!!!!!!!! Allocate arrays !!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(zGnd(soilLevels),porosity(soilLevels),
     +     satPotential(soilLevels),satHydrCond(soilLevels),
     +     soilExponent(soilLevels),heatCapSoil(soilLevels),
     +     u(nsteps),v(nsteps),measRad(nsteps),
     +     scalarFlux(scalarcount))

      allocate(gndScalars(soilLevels,2),scalar(nsteps,scalarcount))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'reading in soil property data'
      
      Open (unit=1,file='inputs/soilTypeParams.ini')
      read(1,*) porosity(1:soilLevels)
      read(1,*) satPotential(1:soilLevels)
      satPotential=satPotential/z_i
      read(1,*) satHydrCond(1:soilLevels)
      satHydrCond=satHydrCond/uScale
      read(1,*) soilExponent(1:soilLevels)
      read(1,*) heatCapSoil(1:soilLevels)
      heatCapSoil=heatCapSoil/(densityAir*Cp_air)
      close(1)
      
      write(*,*) 'reading in LSM parameter info'

      Open (unit=1,file='inputs/soilLevels.ini')
      read(1,*) zGnd(1:soilLevels)
      zGnd=zGnd/z_i
      close(1)

      write(*,*) 'reading in initial soil state'

      Open (unit=1,file='inputs/soilTemperature.ini')
      read(1,*) gndScalars(1:soilLevels,1)
      close(1)
      Open (unit=1,file='inputs/soilMoisture.ini')
      read(1,*) gndScalars(1:soilLevels,2)
      close(1)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'reading in external atmospheric data'
!!  time,u,v,w,theta,q, and tke_sgs at 
!!  z_m for momentum, z_s for scalars 
!!! All of the data are in physical units in the file.   
!!! Also note that time, w, and tkesgs are only read they are not stored.
      
      Open (unit=1,file='inputs/timeseries_10.dat')
      do i=1,nsteps
         read(1,*) time,u(i),v(i),w,scalar(i,1),
     +        scalar(i,2),tkesgs
      enddo
!!!   normalize the input parameters (to match expectation of solveGroundBC)
      u=u/uScale
      v=v/uScale
      do i=1,scalarcount
         scalar(:,i)=scalar(:,i)/scalarScales(i)
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Open files for outputting. Two files, surfaceflux.out for the
!!! terms in the SEB and groundScalar.out for the soil column
!!! profiles

      Open (unit=1,file='outputs/surfaceFlux.out', status='UNKNOWN')
      Open (unit=2,file='outputs/groundScalar.out',status='UNKNOWN')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Call the SEB function for each data point 

      do t = 1,nsteps

!!! Calculate the time (UTC) to use in the radiation model
         UTC = startUTC + float(t)*dt

         M=sqrt(u(t)**2+v(t)**2)

         call solveGroundBC(M,scalar(t,:),gndScalars,ustar,
     +        scalarFlux,soilHeatFlux,porosity,
     +        satPotential,satHydrCond,soilExponent,
     +        heatCapSoil,measRad,
     +        netRad,Psi,Psi0,fi,fi_H,zGnd)
         
!!! write out the surface flux and ground scalar values, note assumes
!!! that we only have 2 scalarflux outputs (sensible and latent). 
!!! This is only for outputting and simply a result of laziness. All
!!! the outputs are dimensional kinematic fluxes and time is in seconds.
!!! to turn the outputs into Watts/m^2 multiply by rho*c_p

         write(1,*) UTC*(z_i/uScale),ustar*uScale,
     +        scalarFlux(1)*scalarScales(1)*uScale,
     +        scalarFlux(2)*scalarScales(2)*uScale,
     +        soilHeatFlux*scalarScales(1)*uScale,
     +	      netRad*scalarScales(1)*uScale

!!! ground scalars are dimensional with units of Kelvin and volumetric content
!!! for temp and moisture, respectively. (Note assumption here that ScalarScales(:)=1
         do i = 1,scalarcount
            write(2,*) (gndScalars(j,i), j=1,soilLevels)
         enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Close the output files
      close(1)
      close(2)

      end
