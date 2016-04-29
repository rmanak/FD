      subroutine ode_driver(lna,lnalpha,a,alpha,Ttt,Trr,Nx,hx,x)
      implicit none
      integer Nx
      real*8 Ttt(Nx),Trr(Nx)    
      real*8 lna(Nx),lnalpha(Nx)
      real*8 x(Nx)
      real*8 ra
      real*8 hx
      real*8 a(Nx),alpha(Nx)
      real*8 lnap1
      real*8 lnalpham1
      real*8 ralpha
      integer i
      integer j,maxiter
      a(1) = 1.0D0
      lna(1) = 0.0D0
      maxiter=10000
      do i=1, Nx-1, 1
      j=0
      lna(i+1) = lna(i)
      ra = 1.0
100   continue
      !/FD generated 1 iteration solver:
      include 'calc_a.inc'
       lna(i+1) = lnap1
      !/FD generated residual:
      include 'resid_a.inc'
       j = j + 1
      if (j .gt. maxiter) then 
         goto 200
      endif
      if (ra .gt. 1.0D-8) then 
          goto 100 
      endif
      a(i+1) = dexp(lna(i+1))
      enddo !/End for
200   continue
      if (j .gt. maxiter) then
          write(0,*) 'iteration didn not converge'
          stop
      end if

      alpha(Nx) = 1.0D0/a(Nx)
      lnalpha(Nx) = dlog(alpha(Nx))
      do i=Nx,2, -1
      j = 0
      lnalpha(i-1) = lnalpha(i)
      ralpha = 1.0D0
300   continue
      include 'calc_alpha.inc'

      lnalpha(i-1) = lnalpham1
      include 'resid_alpha.inc'
      j = j + 1
      if (j .gt. maxiter) then 
         goto 200
      endif
      if (ra .gt. 1.0D-8) then 
          goto 100 
      endif
      alpha(i-1) = dexp(lnalpha(i-1))
      enddo !/End for
400   if (j .gt. maxiter) then
         write(0,*) 'Iteration did not converge'
         stop
      endif

      END

