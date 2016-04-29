      subroutine calc_trr(n_pi,n_pp,Nx,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 n_pi(Nx)
      real*8 n_pp(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      include 'tvd.inc'
      do i=1, Nx, 1
      qb = n_pp(i) ** 2 + n_pi(i) ** 2
      res(i)=qb
      end do
      END
