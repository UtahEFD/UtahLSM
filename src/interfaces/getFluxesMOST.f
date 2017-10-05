      subroutine getFluxesMOST(ustar,Uref,scalarFlux,scalarRef,
     +     Psi,Psi0,fi,PsiH,PsiH0,fiH,computeLH) 
      integer*4 computeLH
      real*8  ustar,Uref,atmTemp,mixingRatio,Psi,Psi0,fi,PsiH,PsiH0,fiH
      real*8, dimension(:) :: scalarFlux, scalarRef
      end subroutine getFluxesMOST
