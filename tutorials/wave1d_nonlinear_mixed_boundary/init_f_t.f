      subroutine init_f_t(v,x,Nx,A,delx,idsignum,x0,res)
      implicit none
      integer i
      integer Nx
      real*8 A
      real*8 delx
      real*8 idsignum
      real*8 x0
      real*8 v(Nx)
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = -0.2D1 * idsignum * v(i) * A * (x(i) - 0.1D1 * x0) / delx ** 
     #2 * exp(-0.1D1 * (x(i) - 0.1D1 * x0) ** 2 / delx ** 2)
      res(i)=qb
      end do
      END
