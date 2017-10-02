      subroutine solveTridiagonalSystem(eDiag,fDiag,gDiag,bDiag,x)
!     given e, f, and g components of tridiagonal matrix (as described below)                                                                                                                           
!     and b (the RHS of a system of equations)                                                                                                                                                          
!     solveTridiagonalSystem() uses the Thomas algorithm to solve the tridiagonal system of eqtns                                                                                                        
!                                                                                                                                                                                                        
!     inputs:                                                                                                                                                                                            
!           e, f, and g to store diagonal 'columns' of data rather than matrix M                                                                                                                         
!           this elimates the unnecessary storage of zeros in M                                                                                                                                          
!           where M would be a tridiagonal matrix, e, f, and g are:                                                                                                                                      
!           matrix is of form, M = [f(1), g(1), 0...       ...0;                                                                                                                                         
!                                  e(2), f(2), g(2), 0... ...0;                                                                                                                                          
!                                  ...                     ...;                                                                                                                                          
!                                  0...        ...0, e(end), f(end)]                                                                                                                                     
!           b - RHS of system of equations where Mx=b, where x are unknowns                                                                                             
!                                                                                                                                                    
!     ouputs:                                                                                                                                             
!           x - array of size soilLevels-1 = size(b), where Mx=b                                                                                                                                         
!               x is the array of unknowns given a tridiagonal system of eqtns                                                                                                                           
      implicit none

      real*8, dimension(:):: bDiag, eDiag, fDiag, gDiag, x

      integer i, length

      length = size(bDiag)

!     decomposition                                                                                                                                                                 
      do i=2,length
         eDiag(i) = eDiag(i)/fDiag(i-1)
         fDiag(i) = fDiag(i) - eDiag(i)*gDiag(i-1)
      enddo

!     forward substitution                                                                                                                                                                     
      do i=2,length
         bDiag(i) = bDiag(i) - eDiag(i)*bDiag(i-1)
      enddo

!     back substitution                                                                                                                                                                              
      x = bDiag/fDiag
      do i = length, 1, -1
         x(i) = ( bDiag(i) - gDiag(i)*x(i+1) )/fDiag(i)
      enddo

      end subroutine solveTridiagonalSystem
