      subroutine ire_f(n_f,nm1_f,np1_f,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 n_f(Nx)
      real*8 nm1_f(Nx)
      real*8 np1_f(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = (nm1_f(i) - 0.2D1 * n_f(i) + np1_f(i)) / ht ** 2 - 0.1D1 * (n
     #_f(i - 1) - 0.2D1 * n_f(i) + n_f(i + 1)) / hx ** 2
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
