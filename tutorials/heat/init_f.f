      subroutine init_f(x,Nx,T0,T1,xmax,xmin,res)
      implicit none
      integer i
      integer Nx
      real*8 T0
      real*8 T1
      real*8 xmax
      real*8 xmin
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = T0 + (T1 - 0.1D1 * T0) * (x(i) - 0.1D1 * xmin) ** 2 / (xmax -
     # 0.1D1 * xmin) ** 2
      res(i)=qb
      end do
      END
