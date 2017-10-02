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
      real*8 minAlbedo,albedo,M,w,tkesgs,time,zo,Psi,
     +     Psi0,fi,fi_H,zlevel,ustar,soilHeatFlux,netRad

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
      read(1,*) satHydrCond(1:soilLevels)
      read(1,*) soilExponent(1:soilLevels)
      read(1,*) heatCapSoil(1:soilLevels)
      close(1)
      
      Open (unit=1,file='inputs/albedo.ini')
      read(1,*) albedo
      close(1)
      
      Open (unit=1,file='inputs/zo.ini')
      read(1,*) zo
      close(1)
      
      write(*,*) 'reading in LSM parameter info'
      
      Open (unit=1,file='inputs/zGnd.ini')
      read(1,*) zGnd(1:soilLevels)
!!! note zGnd is normalized here, all other soil inputs are assumed
!!! normalized before being input (sorry)
      zGnd=zGnd/z_i
      close(1)
      
      if(albedoFlag.eq.1)then
         Open (unit=1,file='inputs/minAlbedo.ini')
         read(1,*) minAlbedo
         close(1)
      endif

      write(*,*) 'reading in initial soil state'

      Open (unit=1,file='inputs/soilTemperature.ini')
      read(1,*) gndScalars(1:soilLevels,1)
      close(1)
      Open (unit=1,file='inputs/soilMoisture.ini')
      read(1,*) gndScalars(1:soilLevels,2)
      close(1)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'reading in external atmospheric data'
!!! Note that currently the file uses timesereis output from LES. 
!!! The time series data file contains ascii output at 10 second intervals
!!! starting at 00:00 UTC of (in order) time,u,v,w,theta,q, and tke_sgs at 
!!! 10 meter hieght (=zlevel) from a single node in the domain (roughly the center).  
!!! All of the data are in physical units in the file.  Note the read in 
!!! could be done in the SEB loop but has been kept separate for clarity.  
!!! Also note that time, w, and tkesgs are only read they are not stored.

      Open (unit=1,file='inputs/zlevel.ini')
      read(1,*) zlevel
      zlevel=zlevel/z_i
      close(1)
      
      Open (unit=1,file='inputs/timeseries_10.dat')
      do i=1,nsteps
         read(1,*) time,u(i),v(i),w,scalar(i,1),
     +        scalar(i,2),tkesgs
      enddo
!!!   normalize the input parameters (to match expectation of solveGroundBC)
      u=u/u_star
      v=v/u_star
      do i=1,scalarcount
         scalar(:,i)=scalar(:,i)/scalarScales(i)
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Open files for outputting. Two files, surfaceflux.out for the
!!! terms in the SEB and groundScalar.out for the soil column
!!! profiles

      Open (unit=1,file='outputs/surfaceflux.out',status='UNKNOWN')
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
     +        heatCapSoil,albedo,minAlbedo,measRad,
     +        netRad,Psi,Psi0,fi,fi_H,zlevel,zo,zGnd)
         
!!! write out the surface flux and ground scalar values, note assumes
!!! that we only have 2 scalarflux outputs (sensible and latent). 
!!! This is only for outputting and simply a result of laziness. All
!!! the outputs are dimensional kinematic fluxes and time is in seconds.
!!! to turn the outputs into Watts/m^2 multiply by rho*c_p

         write(1,*) UTC*(z_i/u_star),ustar*u_star,
     +        scalarFlux(1)*scalarScales(1)*u_star,
     +        scalarFlux(2)*scalarScales(2)*u_star,
     +        soilHeatFlux*scalarScales(1)*u_star,
     +	      netRad*scalarScales(1)*u_star

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
