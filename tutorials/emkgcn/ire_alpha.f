      subroutine ire_alpha(a,Trr,alpha,x,Nx,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 a(Nx)
      real*8 Trr(Nx)
      real*8 alpha(Nx)
      real*8 x(Nx)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = -0.5000000000000000D0 * (alpha(i - 1) - 0.1D1 * alpha(i + 1))
     # / hx / alpha(i) + 0.5000000000000000D0 * (0.1D1 - 0.1D1 * a(i) **
     # 2) / x(i) - 0.5000000000000000D0 * x(i) * Trr(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
