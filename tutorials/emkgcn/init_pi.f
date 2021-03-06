      subroutine init_pi(x,n_phi,Nx,hx,myzero,idsignum,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 idsignum
      real*8 x(Nx)
      real*8 n_phi(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      include 'tvd.inc'
      do i=1, 1, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      do i=2, Nx-1, 1
      qb = idsignum * (0.5000000000000000D0 * (-0.1D1 * n_phi(i - 1) + n
     #_phi(i + 1)) / hx + n_phi(i) / x(i))
      res(i)=qb
      end do
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      END
