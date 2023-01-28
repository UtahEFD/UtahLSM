!!
!! UtahLSM
!! 
!! Copyright (c) 2017–2022 Jeremy A. Gibbs
!! Copyright (c) 2017–2022 Rob Stoll
!! Copyright (c) 2017–2022 Eric Pardyjak
!! Copyright (c) 2017–2022 Pete Willemsen
!! 
!! This file is part of UtahLSM.
!! 
!! This software is free and is distributed under the MIT License.
!! See accompanying LICENSE file or visit https://opensource.org/licenses/MIT.
!!

interface

    function GetInput( input_file ) result( optr )bind(C, name="GetInput")
        import :: c_char, c_ptr
        implicit none
        
        ! Argument list
        character(len=1,kind=c_char), intent(in) :: input_file

        ! Function result
        type(c_ptr) :: optr
    
    end function GetInput
    
    function GetSettings( settings_file ) result( optr )bind(C, name="GetSettings")
        import :: c_char, c_ptr
        implicit none
        
        ! Argument list
        character(len=1,kind=c_char), intent(in) :: settings_file
    
        ! Function result
        type(c_ptr) :: optr
    
    end function GetSettings
    
    function GetOutput( output_file ) result( optr )bind(C, name="GetOutput")
        
        import :: c_char, c_ptr
        implicit none
        
        ! Argument list
        character(len=1,kind=c_char), intent(in) :: output_file
        
        ! Function result
        type(c_ptr) :: optr
    
    end function GetOutput
    
    function GetLSM( input,output,ustar,flux_wT,flux_wq,j,i ) result( optr )bind(C, name="GetLSM")
        
        import :: c_double, c_int, c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: input
        type(c_ptr), value, intent(in) :: output
        real(c_double), intent(in) :: ustar
        real(c_double), intent(in) :: flux_wT
        real(c_double), intent(in) :: flux_wq
        integer(c_int), intent(in) :: j
        integer(c_int), intent(in) :: i
        ! Function result
        type(c_ptr) :: optr
    
    end function GetLSM
    
    subroutine GetItemInt(settings,field,section,name) bind(C, name="GetItemInt")
        
        import :: c_int, c_char, c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: settings
        integer(c_int), intent(in) :: field
        character(len=1,kind=c_char), intent(in) :: section
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetItemInt
    
    subroutine GetDim(input,field,name) bind(C, name="GetDim")
        
        import :: c_int, c_char, c_ptr
        implicit none
         
        ! Argument list
        type(c_ptr), value, intent(in) :: input
        integer(c_int), intent(in) :: field
        character(len=1,kind=c_char), intent(in) :: name
        
    end subroutine GetDim
    
    subroutine GetDataInt(input,field,name) bind(C, name="GetDataInt")
        
        import :: c_int, c_char, c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: input
        integer(c_int), intent(in) :: field
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetDataInt
    
    subroutine GetItemDbl(settings,field,section,name) bind(C, name="GetItemDbl")
        
        import :: c_double, c_char, c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: settings
        real(c_double), intent(in) :: field
        character(len=1,kind=c_char), intent(in) :: section
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetItemDbl
    
    subroutine GetDataDbl(input,field,name) bind(C, name="GetDataDbl")
        
        import :: c_double, c_char, c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: input
        real(c_double), intent(in) :: field
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetDataDbl
    
    subroutine GetItemIntArr(settings,field,size,section,name) bind(C, name="GetItemIntArr")
        
        import :: c_char, c_ptr, c_int
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: settings
        integer(c_int), intent(in) :: field(*)
        integer(c_int), intent(in) :: size
        character(len=1,kind=c_char), intent(in) :: section
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetItemIntArr
    
    subroutine GetDataIntArr(input,field,size,name) bind(C, name="GetDataIntArr")
        
        import :: c_char, c_ptr, c_int
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: input
        integer(c_int), intent(in) :: field(*)
        integer(c_int), intent(in) :: size
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetDataIntArr
    
    subroutine GetItemDblArr(settings,field,size,section,name) bind(C, name="GetItemDblArr")
        
        import :: c_double, c_char, c_ptr, c_int
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: settings
        real(c_double), intent(in) :: field(*)
        integer(c_int), intent(in) :: size
        character(len=1,kind=c_char), intent(in) :: section
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetItemDblArr
    
    subroutine GetDataDblArr(input,field,size,name) bind(C, name="GetDataDblArr")
        
        import :: c_double, c_char, c_ptr, c_int
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: input
        real(c_double), intent(in) :: field(*)
        integer(c_int), intent(in) :: size
        character(len=1,kind=c_char), intent(in) :: name
    
    end subroutine GetDataDblArr
    
    subroutine UpdateFields(lsm,dt,u,T,q,p,rad) bind(C, name="UpdateFields")
        
        import :: c_double, c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: lsm
        real(c_double), intent(in) :: dt
        real(c_double), intent(in) :: u
        real(c_double), intent(in) :: T
        real(c_double), intent(in) :: q
        real(c_double), intent(in) :: p
        real(c_double), intent(in) :: rad
    
    end subroutine UpdateFields
    
    subroutine Run(lsm) bind(C, name="Run")
        
        import :: c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: lsm
    
    end subroutine Run
    
    subroutine Save(lsm,output) bind(C, name="Save")
        
        import :: c_ptr
        implicit none
        
        ! Argument list
        type(c_ptr), value, intent(in) :: lsm
        type(c_ptr), value, intent(in) :: output
    
    end subroutine Save
    
end interface