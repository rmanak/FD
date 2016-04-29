      subroutine res_f(n_f,n_f_t,np1_f,np1_f_t,Nx,ht,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 n_f(Nx)
      real*8 n_f_t(Nx)
      real*8 np1_f(Nx)
      real*8 np1_f_t(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=1, 1, 1
      qb = (-0.1D1 * n_f(i) + np1_f(i)) / ht - 0.5000000000000000D0 * np
     #1_f_t(i) - 0.5000000000000000D0 * n_f_t(i)
      res = res + qb**2
      end do
      do i=2, Nx-1, 1
      qb = (-0.1D1 * n_f(i) + np1_f(i)) / ht - 0.5000000000000000D0 * np
     #1_f_t(i) - 0.5000000000000000D0 * n_f_t(i)
      res = res + qb**2
      end do
      do i=Nx, Nx, 1
      qb = (-0.1D1 * n_f(i) + np1_f(i)) / ht - 0.5000000000000000D0 * np
     #1_f_t(i) - 0.5000000000000000D0 * n_f_t(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
