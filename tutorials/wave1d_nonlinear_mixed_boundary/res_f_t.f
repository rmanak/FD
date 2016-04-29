      subroutine res_f_t(n_f,n_f_t,np1_f,np1_f_t,v,x,Nx,a,b,ht,hx,myzero
     &,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 a
      real*8 b
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 n_f(Nx)
      real*8 n_f_t(Nx)
      real*8 np1_f(Nx)
      real*8 np1_f_t(Nx)
      real*8 v(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=1, 1, 1
      qb = np1_f_t(i) - 0.1D1 * myzero * x(i)
      res = res + qb**2
      end do
      do i=2, Nx-1, 1
      qb = (-0.1D1 * n_f_t(i) + np1_f_t(i)) / ht - 0.5000000000000000D0 
     #* v(i) * (np1_f(i - 1) - 0.2D1 * np1_f(i) + np1_f(i + 1)) / hx ** 
     #2 - 0.5000000000000000D0 * a * np1_f_t(i) - 0.5000000000000000D0 *
     # b * np1_f(i) ** 3 - 0.5000000000000000D0 * v(i) * (n_f(i - 1) - 0
     #.2D1 * n_f(i) + n_f(i + 1)) / hx ** 2 - 0.5000000000000000D0 * a *
     # n_f_t(i) - 0.5000000000000000D0 * b * n_f(i) ** 3
      res = res + qb**2
      end do
      do i=Nx, Nx, 1
      qb = (-0.1D1 * n_f_t(i) + np1_f_t(i)) / ht + 0.2500000000000000D0 
     #* v(i) * (np1_f_t(i - 2) - 0.4D1 * np1_f_t(i - 1) + 0.3D1 * np1_f_
     #t(i)) / hx + 0.2500000000000000D0 * v(i) * (n_f_t(i - 2) - 0.4D1 *
     # n_f_t(i - 1) + 0.3D1 * n_f_t(i)) / hx
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
