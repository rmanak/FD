      subroutine ire_f(np1_f,n_f,nm1_f,Nx,Ny,ht,hx,hy,res)
      implicit none
      integer i
      integer j
      integer Nx
      integer Ny
      real*8 ht
      real*8 hx
      real*8 hy
      real*8 np1_f(Nx,Ny)
      real*8 n_f(Nx,Ny)
      real*8 nm1_f(Nx,Ny)
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=2, Nx-1, 1
      do j=2, Ny-1, 1
      qb = (nm1_f(i, j) - 0.2D1 * n_f(i, j) + np1_f(i, j)) / ht ** 2 + (
     #-0.1D1 * n_f(i - 1, j) + 0.2D1 * n_f(i, j) - 0.1D1 * n_f(i + 1, j)
     #) / hx ** 2 + (-0.1D1 * n_f(i, j - 1) + 0.2D1 * n_f(i, j) - 0.1D1 
     #* n_f(i, j + 1)) / hy ** 2
      res = res + qb**2
      end do
      end do
      res = sqrt(res/(1*Nx*Ny))
      END
