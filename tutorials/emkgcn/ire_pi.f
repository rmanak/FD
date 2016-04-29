      subroutine ire_pi(a,n_pp,nm1_pi,np1_pi,alpha,x,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 a(Nx)
      real*8 n_pp(Nx)
      real*8 nm1_pi(Nx)
      real*8 np1_pi(Nx)
      real*8 alpha(Nx)
      real*8 x(Nx)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = -0.5000000000000000D0 * (nm1_pi(i) - 0.1D1 * np1_pi(i)) / ht 
     #- 0.1D1 / x(i) ** 2 * (0.2D1 * x(i) * alpha(i) / a(i) * n_pp(i) - 
     #0.5000000000000000D0 * x(i) ** 2 * (alpha(i - 1) - 0.1D1 * alpha(i
     # + 1)) / hx / a(i) * n_pp(i) + 0.5000000000000000D0 * x(i) ** 2 * 
     #alpha(i) / a(i) ** 2 * n_pp(i) * (a(i - 1) - 0.1D1 * a(i + 1)) / h
     #x + 0.5000000000000000D0 * x(i) ** 2 * alpha(i) / a(i) * (-0.1D1 *
     # n_pp(i - 1) + n_pp(i + 1)) / hx)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
