!!
!! UtahLSM
!! 
!! Copyright (c) 2017–2023 Jeremy A. Gibbs
!! Copyright (c) 2017–2023 Rob Stoll
!! Copyright (c) 2017–2023 Eric Pardyjak
!! Copyright (c) 2017–2023 Pete Willemsen
!! 
!! This file is part of UtahLSM.
!! 
!! This software is free and is distributed under the MIT License.
!! See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
!! 
program p

   use, intrinsic :: iso_fortran_env, only : compiler_version
   use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_null_ptr, c_associated, c_double, c_int
   
   include "utahlsm_interface.f90"
 
   ! input/output filename variables
   character*50 :: input_offline_file = 'lsm_offline.nc'//char(0)
   character*50 :: settings_file = 'lsm_namelist.json'//char(0)
   character*50 :: input_file = 'lsm_init.nc'//char(0)
   character*50 :: output_file = 'lsm_fortran.nc'//char(0)
   
   ! input/output pointer objects
   type(c_ptr) :: input_offline_obj = c_null_ptr
   type(c_ptr) :: settings_obj = c_null_ptr
   type(c_ptr) :: input_obj = c_null_ptr
   type(c_ptr) :: output_obj = c_null_ptr
   
   ! namelist time section
   character*50 :: name_ntime = 't'//char(0)
   character*50 :: name_tstep = 'tstep'//char(0)
   
   integer :: ntime
   double precision :: tstep
   
   ! namelist time section
   character*50 :: sect_grid = 'grid'//char(0)
   character*50 :: name_imax = 'nx'//char(0)
   character*50 :: name_jmax = 'ny'//char(0)
   
   ! namelist data section
   character*50 :: sect_data = 'data'//char(0)
   character*50 :: name_atmu = 'atm_U'//char(0)
   character*50 :: name_atmt = 'atm_T'//char(0)
   character*50 :: name_atmq = 'atm_q'//char(0)
   character*50 :: name_atmp = 'atm_p'//char(0)
   character*50 :: name_rnet = 'R_net'//char(0)
   
   double precision, allocatable :: atm_U(:),atm_T(:),atm_q(:),atm_p(:),R_net(:)
   
   ! local LSM variables
   double precision :: utc, ustar, flux_wT, flux_wq
   integer :: i, j, k, t
   integer :: i_max, j_max
   type(c_ptr), allocatable :: globalUtahLSM(:)
   
   ! program timing variables
   character*1 :: creturn = achar(13)
   real :: start_time, stop_time, elapsed
   
   ! write a friendly welcome message
   write(*,'(a)')'##############################################################'
   write(*,'(a)')'#                                                            #'
   write(*,'(a)')'#                     Welcome to UtahLSM                     #'
   write(*,'(a)')'#   A land surface model created at the University of Utah   #'
   write(*,'(a)')'#       and the NOAA National Severe Storms Laboratory       #'
   write(*,'(a)')'#                                                            #'
   write(*,'(a)')'##############################################################'
   
   ! Create C++ object representing offline case
   input_offline_obj = GetInput( input_file=input_offline_file )
   
   ! Create C++ object representing model settings
   settings_obj = GetSettings( settings_file=settings_file )
   
   ! Create C++ object representing model init data
   input_obj = GetInput( input_file=input_file )
   
   ! Create C++ object representing model output
   output_obj = GetOutput( output_file=output_file )
      
   ! Get time items from offline input
   call GetDim( input_offline_obj, ntime, name_ntime)
   call GetDataDbl( input_offline_obj, tstep, name_tstep)
   
   ! Get grid items from offline input
   call GetItemInt( settings_obj, i_max, sect_grid, name_imax)
   call GetItemInt( settings_obj, j_max, sect_grid, name_jmax)
   allocate(globalUtahLSM(i_max*j_max))
   
   ! Get items from offline atmospheric data
   allocate(atm_U(ntime))
   allocate(atm_T(ntime))
   allocate(atm_q(ntime))
   allocate(atm_p(ntime))
   allocate(R_net(ntime))
   
   call GetDataDblArr( input_offline_obj, atm_U, ntime, name_atmu)
   call GetDataDblArr( input_offline_obj, atm_T, ntime, name_atmt)
   call GetDataDblArr( input_offline_obj, atm_q, ntime, name_atmq)
   call GetDataDblArr( input_offline_obj, atm_p, ntime, name_atmp)
   call GetDataDblArr( input_offline_obj, R_net, ntime, name_rnet)
   
   ! fill array with LSM instances
   k = 1
   do j = 1, j_max
    do i = 1, i_max
        globalUtahLSM(k) = GetLSM( settings_obj, input_obj, output_obj, ustar, flux_wT, flux_wq, j-1, i-1 )
        k = k+1
    enddo
   enddo
   
   ! set up time information
   call cpu_time(start_time)

   ! loop through all times
   utc = 0
   do t=1,ntime
        utc = utc + tstep
        write(6,'(a,a,f0.2)',advance='no') creturn, '[UtahLSM]        Running for time: ', utc
        flush(6)
    
    ! loop through each LSM instance
    k = 1
    do j = 1, j_max
        do i = 1, i_max
            call UpdateFields(globalUtahLSM(k),tstep,atm_U(t),atm_T(t),atm_q(t),atm_p(t),R_net(t))
            call Run(globalUtahLSM(k))
            call Save(globalUtahLSM(k),output_obj)
            k = k+1
        enddo
    enddo
   enddo
   
   ! compute run time information
   call cpu_time(stop_time)
   elapsed = stop_time - start_time
   write(*,*) creturn
   write(6,'(a,f0.2,a)') "[UtahLSM]        Finished in ",elapsed, " seconds!"
   write(*,'(a)') '##############################################################'

end program p