      subroutine u_f(n_f,x,n_f_t,np1_f,np1_f_t,Nx,Nz,ht,hx,hz,phys_bdy,r
     &es)
      implicit none
      integer i
      integer k
      integer Nx
      integer Nz
      real*8 ht
      real*8 hx
      real*8 hz
      real*8 n_f(Nx,Nz)
      real*8 x(Nx)
      real*8 n_f_t(Nx,Nz)
      real*8 np1_f(Nx,Nz)
      real*8 np1_f_t(Nx,Nz)
      integer phys_bdy(4)
      real*8 res(Nx,Nz)
      real*8 qb
      include 'tvd.inc'
      do i=3, Nx-2, 1
      do k=3, Nz-2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=2, 2, 1
      do k=3, Nz-2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=3, Nz-2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=3, Nx-2, 1
      do k=2, 2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=3, Nx-2, 1
      do k=Nz-1, Nz-1, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=2, 2, 1
      do k=2, 2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=2, 2, 1
      do k=Nz-1, Nz-1, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=2, 2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=Nz-1, Nz-1, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=1, 1, 1
      do k=3, Nz-2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=1, 1, 1
      do k=2, 2, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=1, 1, 1
      do k=Nz-1, Nz-1, 1
      qb = np1_f(i, k) - 0.1D1 * ht * ((-0.1D1 * n_f(i, k) + np1_f(i, k)
     #) / ht - 0.5000000000000000D0 * np1_f_t(i, k) - 0.5000000000000000
     #D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=Nx, Nx, 1
      do k=1, Nz, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hx * x(i) + 0.25D2 * ht * x(
     #i) + 0.6D1 * ht * hx) * ht * hx * x(i) * ((-0.1D1 * n_f(i, k) + np
     #1_f(i, k)) / ht + 0.4166666666666667D-1 * (0.3D1 * np1_f(i - 4, k)
     # - 0.16D2 * np1_f(i - 3, k) + 0.36D2 * np1_f(i - 2, k) - 0.48D2 * 
     #np1_f(i - 1, k) + 0.25D2 * np1_f(i, k)) / hx + 0.2500000000000000D
     #0 * np1_f(i, k) / x(i) + 0.4166666666666667D-1 * (0.3D1 * n_f(i - 
     #4, k) - 0.16D2 * n_f(i - 3, k) + 0.36D2 * n_f(i - 2, k) - 0.48D2 *
     # n_f(i - 1, k) + 0.25D2 * n_f(i, k)) / hx + 0.2500000000000000D0 *
     # n_f(i, k) / x(i))
      res(i,k)=qb
      end do
      end do
      do i=3, Nx-2, 1
      do k=1, 1, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht + 0.4166666666666667D-1
     # * (0.25D2 * np1_f(i, k) - 0.48D2 * np1_f(i, k + 1) + 0.36D2 * np1
     #_f(i, k + 2) - 0.16D2 * np1_f(i, k + 3) + 0.3D1 * np1_f(i, k + 4))
     # / hz - 0.5000000000000000D0 * np1_f_t(i, k) + 0.4166666666666667D
     #-1 * (0.25D2 * n_f(i, k) - 0.48D2 * n_f(i, k + 1) + 0.36D2 * n_f(i
     #, k + 2) - 0.16D2 * n_f(i, k + 3) + 0.3D1 * n_f(i, k + 4)) / hz - 
     #0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=2, 2, 1
      do k=1, 1, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht + 0.4166666666666667D-1
     # * (0.25D2 * np1_f(i, k) - 0.48D2 * np1_f(i, k + 1) + 0.36D2 * np1
     #_f(i, k + 2) - 0.16D2 * np1_f(i, k + 3) + 0.3D1 * np1_f(i, k + 4))
     # / hz - 0.5000000000000000D0 * np1_f_t(i, k) + 0.4166666666666667D
     #-1 * (0.25D2 * n_f(i, k) - 0.48D2 * n_f(i, k + 1) + 0.36D2 * n_f(i
     #, k + 2) - 0.16D2 * n_f(i, k + 3) + 0.3D1 * n_f(i, k + 4)) / hz - 
     #0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=1, 1, 1
      do k=1, 1, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D-1
     # * (-0.25D2 * np1_f(i, k) + 0.48D2 * np1_f(i, k + 1) - 0.36D2 * np
     #1_f(i, k + 2) + 0.16D2 * np1_f(i, k + 3) - 0.3D1 * np1_f(i, k + 4)
     #) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.4166666666666667
     #D-1 * (-0.25D2 * n_f(i, k) + 0.48D2 * n_f(i, k + 1) - 0.36D2 * n_f
     #(i, k + 2) + 0.16D2 * n_f(i, k + 3) - 0.3D1 * n_f(i, k + 4)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=1, 1, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D-1
     # * (-0.25D2 * np1_f(i, k) + 0.48D2 * np1_f(i, k + 1) - 0.36D2 * np
     #1_f(i, k + 2) + 0.16D2 * np1_f(i, k + 3) - 0.3D1 * np1_f(i, k + 4)
     #) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.4166666666666667
     #D-1 * (-0.25D2 * n_f(i, k) + 0.48D2 * n_f(i, k + 1) - 0.36D2 * n_f
     #(i, k + 2) + 0.16D2 * n_f(i, k + 3) - 0.3D1 * n_f(i, k + 4)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=3, Nx-2, 1
      do k=Nz, Nz, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D-1
     # * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2 *
     # np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, k)
     #) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.4166666666666667
     #D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 * 
     #n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=2, 2, 1
      do k=Nz, Nz, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D-1
     # * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2 *
     # np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, k)
     #) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.4166666666666667
     #D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 * 
     #n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=1, 1, 1
      do k=Nz, Nz, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D-1
     # * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2 *
     # np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, k)
     #) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.4166666666666667
     #D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 * 
     #n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=Nz, Nz, 1
      qb = np1_f(i, k) - 0.24D2 / (0.24D2 * hz + 0.25D2 * ht) * ht * hz 
     #* ((-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D-1
     # * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2 *
     # np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, k)
     #) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.4166666666666667
     #D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 * 
     #n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k))
      res(i,k)=qb
      end do
      end do
      END
