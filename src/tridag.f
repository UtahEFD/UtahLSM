      subroutine tridag(a,b,c,r,u,n)
      
      integer*4 n,nmax
      real*8 a(n),b(n),c(n),r(n),u(n)
      
      integer*4 j
      real*8 bet,gam(nmax)
      integer*4, parameter :: nmax=500

      if(b(1).eq.0.) stop 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
           print *,'tridag failed at k=',j
           print *,'a, b, c, gam, and bet=',a(j),b(j),c(j),gam(j),bet     
           stop                   
        end if  
        u(j)=(r(j)-a(j)*u(j-1))/bet
      end do
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      end do
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software