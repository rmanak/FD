      subroutine IRE_a(a,n_phi,nm1_phi,np1_phi,x,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 a(Nx)
      real*8 n_phi(Nx)
      real*8 nm1_phi(Nx)
      real*8 np1_phi(Nx)
      real*8 x(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = -0.5000000000000000D0 * (a(i - 1) - 0.1D1 * a(i + 1)) / hx / 
     #a(i) - 0.5000000000000000D0 * (0.1D1 - 0.1D1 * a(i) ** 2) / x(i) -
     # 0.5000000000000000D0 * x(i) * (0.2500000000000000D0 * (n_phi(i - 
     #1) - 0.1D1 * n_phi(i + 1)) ** 2 / hx ** 2 + 0.2500000000000000D0 *
     # (nm1_phi(i) - 0.1D1 * np1_phi(i)) ** 2 / ht ** 2)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
