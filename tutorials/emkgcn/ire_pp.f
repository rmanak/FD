      subroutine ire_pp(n_pi,nm1_pp,a,alpha,np1_pp,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 n_pi(Nx)
      real*8 nm1_pp(Nx)
      real*8 a(Nx)
      real*8 alpha(Nx)
      real*8 np1_pp(Nx)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_pp(i) + np1_pp(i)) / ht 
     #+ 0.5000000000000000D0 * (alpha(i - 1) - 0.1D1 * alpha(i + 1)) / h
     #x / a(i) * n_pi(i) - 0.5000000000000000D0 * alpha(i) / a(i) ** 2 *
     # n_pi(i) * (a(i - 1) - 0.1D1 * a(i + 1)) / hx - 0.5000000000000000
     #D0 * alpha(i) / a(i) * (-0.1D1 * n_pi(i - 1) + n_pi(i + 1)) / hx
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
