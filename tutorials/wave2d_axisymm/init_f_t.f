      subroutine init_f_t(z,x,Nx,Nz,A,xc,zc,delx,delz,idsigx,idsigz,res)
      implicit none
      integer i
      integer k
      integer Nx
      integer Nz
      real*8 A
      real*8 xc
      real*8 zc
      real*8 delx
      real*8 delz
      real*8 idsigx
      real*8 idsigz
      real*8 z(Nz)
      real*8 x(Nx)
      real*8 res(Nx,Nz)
      real*8 qb
      include 'tvd.inc'
      do i=1, Nx, 1
      do k=1, Nz, 1
      qb = -0.2D1 * idsigx * A * (x(i) - 0.1D1 * xc) / delx ** 2 * exp(-
     #0.1D1 * (x(i) - 0.1D1 * xc) ** 2 / delx ** 2 - 0.1D1 * (z(k) - 0.1
     #D1 * zc) ** 2 / delz ** 2) - 0.2D1 * idsigz * A * (z(k) - 0.1D1 * 
     #zc) / delz ** 2 * exp(-0.1D1 * (x(i) - 0.1D1 * xc) ** 2 / delx ** 
     #2 - 0.1D1 * (z(k) - 0.1D1 * zc) ** 2 / delz ** 2)
      res(i,k)=qb
      end do
      end do
      END
