      subroutine integrateSoilDiffusion(lastSurfScalars,ind)
      use globals
      use SEBmodule
      implicit none

      interface
         include './interfaces/getSoilThermalTransfer.f'
         include './interfaces/getWaterConductivity.f'
         include './interfaces/solveTridiagonalSystem.f'
      end interface

      integer*4 ind
      real*8,dimension(:) :: lastSurfScalars

      integer i
      real*8 dt_
      real*8 D_mid(size(gndScalars,1)-1),D(size(gndScalars,1)),
     >     K(size(gndScalars,1)),dKdz
      real*8 heatCapacity(size(gndScalars,1)),
     >     k_mid(size(gndScalars,1)-1),z_mid(size(gndScalars,1)-1)
!     k_mid and z_mid is the value midway between nodes

!     for moisture diffusion
!     D is diffuse conductivity and K is hydraulic conductivity
      real*8,dimension(size(gndScalars,1)-1) :: b,e,f,g
!     use e, f, and g to store diagonal 'columns' of data rather than a matrix M
!     this elimates the unnecessary storage of zeros in M
!     where M would be a tridiagonal matrix, e, f, and g are:
!     matrix is of form, M = [f(1), g(1), 0...       ...0;
!                             e(2), f(2), g(2), 0... ...0;
!                             ...                     ...;
!                             0...        ...0, e(end), f(end)]


!     actual integration time,
!     not == to dt when soil conditions are not updated every timestep

!     set by integrateSoilDiffFreq
      if( t > endConstSEB )then
         dt_ = dt*integrateSoilDiffFreq
      else
         dt_ = dt
      endif

      if( ind == 1)then
!    compute soil conductivity based on moisture content
         call getSoilThermalTransfer(gndScalars(:,2),K,1,porosity,
     >        satPotential,soilExponent,heatCapSoil)

!     interpolate zGnd and k to get values at mid-levels
!     for z do this once and save z_mid (never changes)
         do i = 1,soilLevels-1
            k_mid(i) = sum(k(i:i+1))/2.d0
            z_mid(i) = sum(zGnd(i:i+1))/2.d0
         enddo

      else

!     compute diffusive and hydraulic water conductivity of soil
!     given the soil properties from getGroundParams.f
         call getWaterConductivity(gndScalars(:,2),D,K,porosity,
     >        satPotential,satHydrCond,soilExponent)

         do i = 1,soilLevels-1
            z_mid(i) = sum(zGnd(i:i+1))/2.d0
            k_mid(i) = sum(D(i:i+1))/2.d0
         enddo
         
      endif

!     set coefficients for node 2, the first row of M (top node, surface is solved )
!     node 2 and soilLevels coefficients have a different form of implicit finite diff.
!     because of the boundary conditions.
      f(1)   = ( dt_/(2.d0*(z_mid(2)-z_mid(1))) ) *
     >     ( k_mid(1)/(zGnd(2)-zGnd(1)) + k_mid(2)/(zGnd(3)-zGnd(2)))+
     >     1.d0
      g(1)   = - dt_*k_mid(2)/(2.d0*(z_mid(2)-z_mid(1))*(zGnd(3)-
     >     zGnd(2)))

      b(1)   = (dt_/(2.d0*(z_mid(2)-z_mid(1))))
     >     * ( (k_mid(2)*(gndScalars(3,ind) - gndScalars(2,ind))
     >     / (zGnd(3)-zGnd(2))) - ( k_mid(1)
     >     * (gndScalars(2,ind) - lastSurfScalars(ind))
     >     /(zGnd(2)-zGnd(1))) )
     >     + gndScalars(2,ind) + gndScalars(1,ind)*dt_*k_mid(1)
     >     / (2.d0*(z_mid(2)-z_mid(1))*(zGnd(2)-zGnd(1)))

      if( ind == 2 )then
         dKdz = dt_*(K(1)*(zGnd(2)-zGnd(3))
     >        / ((zGnd(1)-zGnd(2))*(zGnd(1)-zGnd(3)))
     >        + K(2)*(2.d0*zGnd(2)-zGnd(1)-zGnd(3))
     >        / ((zGnd(2)-zGnd(1))*(zGnd(2)-zGnd(3)))
     >        + K(3)*(zGnd(2)-zGnd(1))
     >        / ((zGnd(3)-zGnd(1))*(zGnd(3)-zGnd(2))))
         
         b(1) = b(1) + dKdz
      endif
      
!     set coefficients for node soilLevels (last ground node)
      e(soilLevels-1) = -dt_*k_mid(soilLevels-1)
     >     / (2.d0*(zGnd(soilLevels)-zGnd(soilLevels-1))**2)
      f(soilLevels-1) = 1.d0 + dt_*k_mid(soilLevels-1)
     >     /(2.d0*(zGnd(soilLevels)-zGnd(soilLevels-1))**2)
      b(soilLevels-1) = gndScalars(soilLevels,ind)
     >     - ( gndScalars(soilLevels,ind)
     >     - gndScalars(soilLevels-1,ind) )*dt_*k_mid(soilLevels-1)
     >     /(2.d0*(zGnd(soilLevels)-zGnd(soilLevels-1))**2)
      
      if( ind == 2 )then
         dKdz = dt_*(K(soilLevels)-K(soilLevels-1))
     >        / (zGnd(soilLevels)-zGnd(soilLevels-1))
         
         b(soilLevels-1) = b(soilLevels-1) + dKdz
      endif

!     set matrix coefficients for nodes 3 through soilLevels-1
!     form of these coefficients is identical, there are no boundary conditions
      if( soilLevels > 3)then
         do i = 3,soilLevels-1
            e(i-1) = -dt_*k_mid(i-1)/(2.d0*(z_mid(i)-z_mid(i-1))*
     >           (zGnd(i)-zGnd(i-1))) !M(i-1,i-2)

            f(i-1) = 1.d0 + dt_*k_mid(i-1)
     >           / (2.d0*(z_mid(i)-z_mid(i-1))*(zGnd(i)-zGnd(i-1)))
     >           +dt_*k_mid(i)
     >           / (2.d0*(z_mid(i)-z_mid(i-1))*(zGnd(i+1)-zGnd(i))) !M(i-1,i-1)

            g(i-1)   = -dt_*k_mid(i) /
     >           (2.d0*(z_mid(i)-z_mid(i-1))*(zGnd(i+1)-zGnd(i))) !M(i-1,i)

            b(i-1)     = gndScalars(i,ind)
     >           + (dt_/(2.d0*(z_mid(i)-z_mid(i-1))))*
     >           (k_mid(i)*(gndScalars(i+1,ind)
     >           -gndScalars(i,ind))/(zGnd(i+1)-zGnd(i))
     >           -k_mid(i-1)*(gndScalars(i,ind)-gndScalars(i-1,ind))
     >           / (zGnd(i)-zGnd(i-1)))

            if( ind == 2)then
!     explicit finite difference for unevenly spaced data
               dKdz = dt_*(K(i-1)*(zGnd(i)-zGnd(i+1))
     >              / ((zGnd(i-1)-zGnd(i))*(zGnd(i-1)-zGnd(i+1)))
     >              + K(i)*(2.d0*zGnd(i)-zGnd(i-1)-zGnd(i+1))
     >              / ((zGnd(i)-zGnd(i-1))*(zGnd(i)-zGnd(i+1)))
     >              + K(i+1)*(zGnd(i)-zGnd(i-1))
     >              / ((zGnd(i+1)-zGnd(i-1))*(zGnd(i+1)-zGnd(i))))
               
               b(i-1) = b(i-1) + dKdz
               
            endif
            
         enddo
      endif
      
c     call solveTridiagonalSystem(e,f,g,b,
c     >     gndScalars(2:size(gndScalars,1),ind))
      call tridag(e,f,g,b,gndScalars(2:soilLevels,ind),soilLevels-1)
      
      end subroutine integrateSoilDiffusion
