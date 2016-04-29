      subroutine init_phi(x,Nx,A,x0,delx,res)
      implicit none
      integer i
      integer Nx
      real*8 A
      real*8 x0
      real*8 delx
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      include 'tvd.inc'
      do i=1, Nx, 1
      qb = A * exp(-0.1D1 * (x(i) - 0.1D1 * x0) ** 2 / delx ** 2)
      res(i)=qb
      end do
      END
