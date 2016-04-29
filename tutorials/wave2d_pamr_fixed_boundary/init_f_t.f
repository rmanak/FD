      subroutine init_f_t(y,x,Nx,Ny,A,xc,yc,delx,dely,idsigx,idsigy,res)
      implicit none
      integer i
      integer j
      integer Nx
      integer Ny
      real*8 A
      real*8 xc
      real*8 yc
      real*8 delx
      real*8 dely
      real*8 idsigx
      real*8 idsigy
      real*8 y(Ny)
      real*8 x(Nx)
      real*8 res(Nx,Ny)
      real*8 qb
      include 'tvd.inc'
      do i=1, Nx, 1
      do j=1, Ny, 1
      qb = -0.2D1 * idsigx * A * (x(i) - 0.1D1 * xc) / delx ** 2 * exp(-
     #0.1D1 * (x(i) - 0.1D1 * xc) ** 2 / delx ** 2 - 0.1D1 * (y(j) - 0.1
     #D1 * yc) ** 2 / dely ** 2) - 0.2D1 * idsigy * A * (y(j) - 0.1D1 * 
     #yc) / dely ** 2 * exp(-0.1D1 * (x(i) - 0.1D1 * xc) ** 2 / delx ** 
     #2 - 0.1D1 * (y(j) - 0.1D1 * yc) ** 2 / dely ** 2)
      res(i,j)=qb
      end do
      end do
      END
