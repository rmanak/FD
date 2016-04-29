      subroutine u_f(n_f,x,y,n_f_t,np1_f,np1_f_t,Nx,Ny,ht,myzero,phys_bd
     &y,res)
      implicit none
      integer i
      integer j
      integer Nx
      integer Ny
      real*8 ht
      real*8 myzero
      real*8 n_f(Nx,Ny)
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 n_f_t(Nx,Ny)
      real*8 np1_f(Nx,Ny)
      real*8 np1_f_t(Nx,Ny)
      integer phys_bdy(4)
      real*8 res(Nx,Ny)
      real*8 qb
      include 'tvd.inc'
      do i=2, Nx-1, 1
      do j=2, Ny-1, 1
      qb = np1_f(i, j) - 0.1D1 * ht * ((-0.1D1 * n_f(i, j) + np1_f(i, j)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, j) - 0.5000000000000000
     #D0 * n_f_t(i, j))
      res(i,j)=qb
      end do
      end do
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, Ny, 1
      qb = myzero * x(i) * y(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      do j=1, Ny, 1
      qb = myzero * x(i) * y(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=1, Nx, 1
      do j=1, 1, 1
      qb = myzero * x(i) * y(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=1, Nx, 1
      do j=Ny, Ny, 1
      qb = myzero * x(i) * y(j)
      res(i,j)=qb
      end do
      end do
      endif
      END
