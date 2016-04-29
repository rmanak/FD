      subroutine eval_laplace(phi,x,z,Nx,Nz,hx,hz,myzero,phys_bdy,res)
      implicit none
      integer i
      integer k
      integer Nx
      integer Nz
      real*8 hx
      real*8 hz
      real*8 myzero
      real*8 phi(Nx,Nz)
      real*8 x(Nx)
      real*8 z(Nz)
      integer phys_bdy(4)
      real*8 res(Nx,Nz)
      real*8 qb
      do i=2, Nx-1, 1
      do k=2, Nz-1, 1
      qb = (phi(i - 1, k) - 0.2D1 * phi(i, k) + phi(i + 1, k)) / hx ** 2
     # - 0.5000000000000000D0 * (phi(i - 1, k) - 0.1D1 * phi(i + 1, k)) 
     #/ hx / x(i) + (phi(i, k - 1) - 0.2D1 * phi(i, k) + phi(i, k + 1)) 
     #/ hz ** 2
      res(i,k)=qb
      end do
      end do
      do i=1, 1, 1
      do k=2, Nz-1, 1
      qb = 0.2D1 * (0.2D1 * phi(i + 1, k) - 0.2D1 * phi(i, k)) / hx ** 2
     # + (phi(i, k - 1) - 0.2D1 * phi(i, k) + phi(i, k + 1)) / hz ** 2
      res(i,k)=qb
      end do
      end do
      do i=Nx, Nx, 1
      do k=1, Nz, 1
      qb = myzero * x(i) * z(k)
      res(i,k)=qb
      end do
      end do
      do i=1, Nx, 1
      do k=1, 1, 1
      qb = myzero * x(i) * z(k)
      res(i,k)=qb
      end do
      end do
      do i=1, Nx, 1
      do k=Nz, Nz, 1
      qb = myzero * x(i) * z(k)
      res(i,k)=qb
      end do
      end do
      END
