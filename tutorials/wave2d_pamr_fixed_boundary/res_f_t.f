      subroutine res_f_t(n_f,x,y,n_f_t,np1_f,np1_f_t,Nx,Ny,ht,hx,hy,myze
     &ro,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer Nx
      integer Ny
      real*8 ht
      real*8 hx
      real*8 hy
      real*8 myzero
      real*8 n_f(Nx,Ny)
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 n_f_t(Nx,Ny)
      real*8 np1_f(Nx,Ny)
      real*8 np1_f_t(Nx,Ny)
      integer phys_bdy(4)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=2, Nx-1, 1
      do j=2, Ny-1, 1
      qb = (-0.1D1 * n_f_t(i, j) + np1_f_t(i, j)) / ht - 0.5000000000000
     #000D0 * (np1_f(i - 1, j) - 0.2D1 * np1_f(i, j) + np1_f(i + 1, j)) 
     #/ hx ** 2 - 0.5000000000000000D0 * (np1_f(i, j - 1) - 0.2D1 * np1_
     #f(i, j) + np1_f(i, j + 1)) / hy ** 2 - 0.5000000000000000D0 * (n_f
     #(i - 1, j) - 0.2D1 * n_f(i, j) + n_f(i + 1, j)) / hx ** 2 - 0.5000
     #000000000000D0 * (n_f(i, j - 1) - 0.2D1 * n_f(i, j) + n_f(i, j + 1
     #)) / hy ** 2
      res = res + qb**2
      end do
      end do
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, Ny, 1
      qb = np1_f_t(i, j) - 0.1D1 * myzero * x(i) * y(j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      do j=1, Ny, 1
      qb = np1_f_t(i, j) - 0.1D1 * myzero * x(i) * y(j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=1, Nx, 1
      do j=1, 1, 1
      qb = np1_f_t(i, j) - 0.1D1 * myzero * x(i) * y(j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=1, Nx, 1
      do j=Ny, Ny, 1
      qb = np1_f_t(i, j) - 0.1D1 * myzero * x(i) * y(j)
      res = res + qb**2
      end do
      end do
      endif
      res = sqrt(res/(1*Nx*Ny))
      END
