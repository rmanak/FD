      subroutine update_f(n_f,np1_f,x,Nx,T0,T1,ht,hx,myzero,phys_bdy,res
     &)
      implicit none
      integer i
      integer Nx
      real*8 T0
      real*8 T1
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 n_f(Nx)
      real*8 np1_f(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      do i=1, 1, 1
      qb = T0 - 0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      do i=2, Nx-1, 1
      qb = np1_f(i) - 0.1D1 * ht * ((-0.1D1 * n_f(i) + np1_f(i)) / ht - 
     #0.1D1 * (n_f(i - 1) - 0.2D1 * n_f(i) + n_f(i + 1)) / hx ** 2)
      res(i)=qb
      end do
      do i=Nx, Nx, 1
      qb = T1 - 0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      END
