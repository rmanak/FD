      subroutine ire_f(n_f,nm2_f,np1_f,nm1_f,x,np2_f,Nx,Nz,ht,hx,hz,res)
      implicit none
      integer i
      integer k
      integer Nx
      integer Nz
      real*8 ht
      real*8 hx
      real*8 hz
      real*8 n_f(Nx,Nz)
      real*8 nm2_f(Nx,Nz)
      real*8 np1_f(Nx,Nz)
      real*8 nm1_f(Nx,Nz)
      real*8 x(Nx)
      real*8 np2_f(Nx,Nz)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=3, Nx-2, 1
      do k=3, Nz-2, 1
      qb = -0.8333333333333333D-1 * (nm2_f(i, k) - 0.16D2 * nm1_f(i, k) 
     #+ 0.30D2 * n_f(i, k) - 0.16D2 * np1_f(i, k) + np2_f(i, k)) / ht **
     # 2 - 0.1D1 / x(i) * (-0.8333333333333333D-1 * (-0.1D1 * n_f(i - 2,
     # k) + 0.8D1 * n_f(i - 1, k) - 0.8D1 * n_f(i + 1, k) + n_f(i + 2, k
     #)) / hx - 0.8333333333333333D-1 * x(i) * (n_f(i - 2, k) - 0.16D2 *
     # n_f(i - 1, k) + 0.30D2 * n_f(i, k) - 0.16D2 * n_f(i + 1, k) + n_f
     #(i + 2, k)) / hx ** 2) - 0.8333333333333333D-1 * (-0.1D1 * n_f(i, 
     #k - 2) + 0.16D2 * n_f(i, k - 1) - 0.30D2 * n_f(i, k) + 0.16D2 * n_
     #f(i, k + 1) - 0.1D1 * n_f(i, k + 2)) / hz ** 2
      res = res + qb**2
      end do
      end do
      res = sqrt(res/(1*Nx*Nz))
      END
