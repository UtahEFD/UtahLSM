!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         This program is developed to run the LSM          !
!         model developed by Shingleton (2010)              !
!                                                           !
!         v1.0: developed by Shingleton and Stoll 2012      !
!         v2.0: developed by Gibbs and Stoll 2017           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Program utahLSM

      use globals
      use SEBmodule
      implicit none
      
      interface
         include './interfaces/solveGroundBC.f'
      end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Variable definition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer*4 i,j
      real*8 M,soilHeatFlux,netRad,ustar
      real*8,allocatable,dimension(:) :: S, scalarFlux
      character(len=32) :: caseName

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Parse command line for input case name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (command_argument_count().ne.1) then
        write(*,*) "Usage: ./LSM case_name"
        write(*,*) "where case_name is subfolder in inputs/"
        stop
      end if
	  
	  call get_command_argument(1, caseName)
	  	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read inputs from LSMinputs and external soil/atm data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'reading inputs'
      call readInputs(caseName)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Allocate arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(S(scalarcount), scalarFlux(scalarcount))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Open output files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Open (unit=1,file='outputs/surfaceFlux.out', status='UNKNOWN')
      Open (unit=2,file='outputs/groundScalar.out',status='UNKNOWN')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Solve SEB in time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do t = 1,nsteps

         ! calculate time (UTC) to use in the radiation model
         UTC = startUTC + float(t)*dt
         
         ! compute magnitude of horizontal wind vector
         M=sqrt(u(t)**2+v(t)**2)
         
         ! grab atmospheric scalars
         S=scalar(t,:)
         
         ! call the model
         call solveGroundBC(M,S,ustar,scalarFlux,soilHeatFlux,netRad)
     
         ! write fluxes in dimensional kinematic units
         write(1,*) UTC*(z_i/uScale),ustar*uScale,
     +        scalarFlux(1)*scalarScales(1)*uScale,
     +        scalarFlux(2)*scalarScales(2)*uScale,
     +        soilHeatFlux*scalarScales(1)*uScale,
     +        netRad*scalarScales(1)*uScale

        ! write soil scalars in units of K and volum. content
         do i = 1,scalarcount
            write(2,*) (gndScalars(j,i), j=1,soilLevels)
         enddo

      enddo

      ! Close the output files
      close(1)
      close(2)

      end