      subroutine u_f_t(n_f,n_f_t,np1_f,np1_f_t,Nx,ht,hx,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 n_f(Nx)
      real*8 n_f_t(Nx)
      real*8 np1_f(Nx)
      real*8 np1_f_t(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      do i=1, 1, 1
      qb = np1_f_t(i) - 0.1D1 * ht * ((-0.1D1 * n_f_t(i) + np1_f_t(i)) /
     # ht - 0.5000000000000000D0 * (np1_f(i - 2 + Nx) - 0.2D1 * np1_f(i)
     # + np1_f(i + 1)) / hx ** 2 - 0.5000000000000000D0 * (n_f(i - 2 + N
     #x) - 0.2D1 * n_f(i) + n_f(i + 1)) / hx ** 2)
      res(i)=qb
      end do
      do i=2, Nx-1, 1
      qb = np1_f_t(i) - 0.1D1 * ht * ((-0.1D1 * n_f_t(i) + np1_f_t(i)) /
     # ht - 0.5000000000000000D0 * (np1_f(i - 1) - 0.2D1 * np1_f(i) + np
     #1_f(i + 1)) / hx ** 2 - 0.5000000000000000D0 * (n_f(i - 1) - 0.2D1
     # * n_f(i) + n_f(i + 1)) / hx ** 2)
      res(i)=qb
      end do
      do i=Nx, Nx, 1
      qb = np1_f_t(i) - 0.1D1 * ht * ((-0.1D1 * n_f_t(i) + np1_f_t(i)) /
     # ht - 0.5000000000000000D0 * (np1_f(i - 1) - 0.2D1 * np1_f(i) + np
     #1_f(i + 2 - Nx)) / hx ** 2 - 0.5000000000000000D0 * (n_f(i - 1) - 
     #0.2D1 * n_f(i) + n_f(i + 2 - Nx)) / hx ** 2)
      res(i)=qb
      end do
      END
