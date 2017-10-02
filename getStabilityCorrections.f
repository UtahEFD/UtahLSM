      subroutine getStabilityCorrections(zlevel,ustar,atmTemp,
     +     mixingRatio,scalarFlux,Psi,Psi0,fi,PsiH,PsiH0,fiH)
      
      use globals
      use SEBmodule
      implicit none
      
      real*8 zlevel,ustar,atmTemp,mixingRatio,Psi,Psi0,fi,PsiH,PsiH0,fiH
      real*8, dimension(:) :: scalarFlux
      
      real*8 tempFlux,obukhovLength,x,y
      
      tempFlux = atmTemp*0.61d0*scalarFlux(2)
     >     + scalarFlux(1)*(1.d0+0.61d0*mixingRatio)
      
      obukhovLength = -ustar**3*atmTemp /
     >     ( vonk*g_hat*tempFlux )
      
!!! the cap on stability can probably be lifted here when not LES !!!                                                                                                               
      if ( zlevel/obukhovLength .gt. 5.0 )then
         obukhovLength = zlevel/5.0
      elseif( zlevel/obukhovLength .lt. -5.0)then
         obukhovLength = -zlevel/5.0
      endif
      
      if( tempFlux.gt.0.0 )then
         
         x=(1.d0-(15.d0*zlevel/obukhovLength))**0.25d0
         y=(1.d0-(15.d0*(zt/z_i)/obukhovLength))**0.25d0 ! zo changed to zt           
         Psi = -2.d0*dlog(0.5d0*(1+x)) -
     >        dlog(0.5d0*(1.d0+x**2.d0))+2.d0*atan(x)-pi/2.d0
         Psi0 = -2.d0*dlog(0.5d0*(1+y)) -
     >        dlog(0.5d0*(1.d0+y**2.d0)) + 2.d0*atan(y) - pi/2.d0
         fi = 1.d0/x
!     eqtn 11.9 and 11.14 from Arya                                                                                                                                                   
         PsiH  = -2.d0*dlog(0.5d0*(1.d0+x**2.d0))
         PsiH0 = -2.d0*dlog(0.5d0*(1+y**2.d0))
         fiH   = fi**2.d0
         
      elseif( tempFlux.lt.-0.0 )then
         
!     Psi   = 4.7*z/obukhovLength                                                                                                                                                                
!     for GABLS 5*z/L usually 4.7*z/L                                                                                                                                       
         
         Psi   = 5.d0*zlevel/obukhovLength
         Psi0  = 5.d0*(zt/z_i)/obukhovLength !zo changed to zt                                                                                                                  
         fi    = 1.d0 + Psi
         PsiH  = 5.d0*zlevel/obukhovLength
         PsiH0 = 5.d0*(zt/z_i)/obukhovLength ! zo changed to zt                                      
         fiH   = 0.74d0 + PsiH
      else
         
         Psi   = 0.d0
         Psi0  = 0.d0
         fi    = 1.d0
         PsiH  = Psi
         PsiH0 = Psi0
         fiH   = 0.74d0
         
      endif
      
      end subroutine getStabilityCorrections
      
