      subroutine init_v(x,Nx,d,p,vcons,res)
      implicit none
      integer i
      integer Nx
      real*8 d
      real*8 p
      real*8 vcons
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.10D1 + vcons * tanh((x(i) - 0.1D1 * p) / d)
      res(i)=qb
      end do
      END
