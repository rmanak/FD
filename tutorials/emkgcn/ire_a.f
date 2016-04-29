      subroutine ire_a(Ttt,a,x,Nx,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 Ttt(Nx)
      real*8 a(Nx)
      real*8 x(Nx)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = -0.5000000000000000D0 * (a(i - 1) - 0.1D1 * a(i + 1)) / hx / 
     #a(i) - 0.5000000000000000D0 * (0.1D1 - 0.1D1 * a(i) ** 2) / x(i) +
     # 0.5000000000000000D0 * x(i) * Ttt(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
