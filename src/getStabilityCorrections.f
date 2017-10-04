      subroutine getStabilityCorrections(ustar,atmTemp,
     +     mixingRatio,scalarFlux,Psi,Psi0,fi,PsiH,PsiH0,fiH)
      
      use globals
      use SEBmodule
      implicit none
      
      real*8  ustar,atmTemp,mixingRatio,Psi,Psi0,fi,PsiH,PsiH0,fiH
      real*8, dimension(:) :: scalarFlux
      
      real*8 tempFlux,obukhovLength,x,y
      
      tempFlux = atmTemp*0.61d0*scalarFlux(2)
     >     + scalarFlux(1)*(1.d0+0.61d0*mixingRatio)
      
      obukhovLength = -ustar**3*atmTemp /
     >     ( vonk*g_hat*tempFlux )
      
      ! the cap can probably be lifted here when not LES !!!                                                                                                               
      !if ( z_m/obukhovLength .gt. 5.0 )then
      !   obukhovLength = z_m/5.0
      !elseif( z_m/obukhovLength .lt. -5.0)then
      !   obukhovLength = -z_m/5.0
      !endif
     
      ! unstable 
      if( tempFlux.gt.0.0 )then
         
         ! momentum         
         x=(1.d0-(15.d0*z_m/obukhovLength))**0.25d0
         y=(1.d0-(15.d0*z_o/obukhovLength))**0.25d0           
         Psi = -2.d0*dlog(0.5d0*(1+x)) -
     >        dlog(0.5d0*(1.d0+x**2.d0))+2.d0*atan(x)-pi/2.d0
         Psi0 = -2.d0*dlog(0.5d0*(1+y)) -
     >        dlog(0.5d0*(1.d0+y**2.d0)) + 2.d0*atan(y) - pi/2.d0
         fi = 1.d0/x

         ! scalar    
         x=(1.d0-(15.d0*z_s/obukhovLength))**0.25d0
         y=(1.d0-(15.d0*z_t/obukhovLength))**0.25d0                                                                                                                                              
         PsiH  = -2.d0*dlog(0.5d0*(1.d0+x**2.d0))
         PsiH0 = -2.d0*dlog(0.5d0*(1.d0+y**2.d0))
         fiH   = fi**2.d0

      ! stable
      elseif( tempFlux.lt.-0.0 )then
         
         ! momentum                                                                                                                                     
         Psi   = 5.d0*z_m/obukhovLength
         Psi0  = 5.d0*z_o/obukhovLength                                                                                                                
         fi    = 1.d0 + Psi
         
         ! scalar
         PsiH  = 5.d0*z_s/obukhovLength
         PsiH0 = 5.d0*z_t/obukhovLength                                    
         fiH   = 0.74d0 + PsiH
      
      ! neutral
      else
         
         ! momentum
         Psi   = 0.d0
         Psi0  = 0.d0
         fi    = 1.d0
         
         ! scalar
         PsiH  = Psi
         PsiH0 = Psi0
         fiH   = 0.74d0
         
      endif
      
      end subroutine getStabilityCorrections