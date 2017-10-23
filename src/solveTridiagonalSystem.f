      subroutine solveTridiagonalSystem(eDiag,fDiag,gDiag,bDiag,x)
      ! This function uses the Thomas algorithm to solve the tridiagonal
      ! system of equations                                                                                                        
      !                                                                                                                                                                                                    
      ! inputs:                                                                                                                                                                                            
      !       e, f, and g store diagonal 'columns' of data rather than 
      !       matrix M to elimate the unnecessary storage of zeros in M,
      !       where M is a tridiagonal matrix of the form:                                                                                                                                      
      !           M = [f(1), g(1), 0...       ...0;                                                                                                                                         
      !                e(2), f(2), g(2), 0... ...0;                                                                                                                                          
      !                ...                     ...;                                                                                                                                          
      !                0...        ...0, e(end), f(end)]                                                                                                                                     
      !       b - RHS of system of equations Mx=b, where x are unknowns
      !                                                                                                                                                
      ! ouputs:                                                                                                                                             
      !       x - array of size soilLevels-1 = size(b), where Mx=b and                                                                                                                                         
      !           x is the array of unknowns given a tridiagonal system
      !           of equations                                                                                                                           
      implicit none

      real*8, dimension(:):: bDiag, eDiag, fDiag, gDiag, x

      integer i, length

      length = size(bDiag)

      ! decomposition                                                                                                                                                                 
      do i=2,length
         eDiag(i) = eDiag(i)/fDiag(i-1)
         fDiag(i) = fDiag(i) - eDiag(i)*gDiag(i-1)
      enddo

      ! forward substitution                                                                                                                                                                     
      do i=2,length
         bDiag(i) = bDiag(i) - eDiag(i)*bDiag(i-1)
      enddo

      ! back substitution                                                                                                                                                                              
      x = bDiag/fDiag
      do i = length, 1, -1
         x(i) = ( bDiag(i) - gDiag(i)*x(i+1) )/fDiag(i)
      enddo

      end subroutine solveTridiagonalSystem
