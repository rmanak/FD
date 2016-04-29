      subroutine res_f(n_f,x,n_f_t,np1_f,np1_f_t,Nx,Nz,ht,hx,hz,phys_bdy
     &,res)
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
      real*8 res
      real*8 qb
      include 'tvd.inc'
      res = 0.0D0
      do i=3, Nx-2, 1
      do k=3, Nz-2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=2, 2, 1
      do k=3, Nz-2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=3, Nz-2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=3, Nx-2, 1
      do k=2, 2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=3, Nx-2, 1
      do k=Nz-1, Nz-1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=2, 2, 1
      do k=2, 2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=2, 2, 1
      do k=Nz-1, Nz-1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=2, 2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=Nz-1, Nz-1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=1, 1, 1
      do k=3, Nz-2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=1, 1, 1
      do k=2, 2, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=1, 1, 1
      do k=Nz-1, Nz-1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, k) - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=Nx, Nx, 1
      do k=1, Nz, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht + 0.4166666666666667D
     #-1 * (0.3D1 * np1_f(i - 4, k) - 0.16D2 * np1_f(i - 3, k) + 0.36D2 
     #* np1_f(i - 2, k) - 0.48D2 * np1_f(i - 1, k) + 0.25D2 * np1_f(i, k
     #)) / hx + 0.2500000000000000D0 * np1_f(i, k) / x(i) + 0.4166666666
     #666667D-1 * (0.3D1 * n_f(i - 4, k) - 0.16D2 * n_f(i - 3, k) + 0.36
     #D2 * n_f(i - 2, k) - 0.48D2 * n_f(i - 1, k) + 0.25D2 * n_f(i, k)) 
     #/ hx + 0.2500000000000000D0 * n_f(i, k) / x(i)
      res = res + qb**2
      end do
      end do
      do i=3, Nx-2, 1
      do k=1, 1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht + 0.4166666666666667D
     #-1 * (0.25D2 * np1_f(i, k) - 0.48D2 * np1_f(i, k + 1) + 0.36D2 * n
     #p1_f(i, k + 2) - 0.16D2 * np1_f(i, k + 3) + 0.3D1 * np1_f(i, k + 4
     #)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) + 0.416666666666666
     #7D-1 * (0.25D2 * n_f(i, k) - 0.48D2 * n_f(i, k + 1) + 0.36D2 * n_f
     #(i, k + 2) - 0.16D2 * n_f(i, k + 3) + 0.3D1 * n_f(i, k + 4)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=2, 2, 1
      do k=1, 1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht + 0.4166666666666667D
     #-1 * (0.25D2 * np1_f(i, k) - 0.48D2 * np1_f(i, k + 1) + 0.36D2 * n
     #p1_f(i, k + 2) - 0.16D2 * np1_f(i, k + 3) + 0.3D1 * np1_f(i, k + 4
     #)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) + 0.416666666666666
     #7D-1 * (0.25D2 * n_f(i, k) - 0.48D2 * n_f(i, k + 1) + 0.36D2 * n_f
     #(i, k + 2) - 0.16D2 * n_f(i, k + 3) + 0.3D1 * n_f(i, k + 4)) / hz 
     #- 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=1, 1, 1
      do k=1, 1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D
     #-1 * (-0.25D2 * np1_f(i, k) + 0.48D2 * np1_f(i, k + 1) - 0.36D2 * 
     #np1_f(i, k + 2) + 0.16D2 * np1_f(i, k + 3) - 0.3D1 * np1_f(i, k + 
     #4)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.41666666666666
     #67D-1 * (-0.25D2 * n_f(i, k) + 0.48D2 * n_f(i, k + 1) - 0.36D2 * n
     #_f(i, k + 2) + 0.16D2 * n_f(i, k + 3) - 0.3D1 * n_f(i, k + 4)) / h
     #z - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=1, 1, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D
     #-1 * (-0.25D2 * np1_f(i, k) + 0.48D2 * np1_f(i, k + 1) - 0.36D2 * 
     #np1_f(i, k + 2) + 0.16D2 * np1_f(i, k + 3) - 0.3D1 * np1_f(i, k + 
     #4)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.41666666666666
     #67D-1 * (-0.25D2 * n_f(i, k) + 0.48D2 * n_f(i, k + 1) - 0.36D2 * n
     #_f(i, k + 2) + 0.16D2 * n_f(i, k + 3) - 0.3D1 * n_f(i, k + 4)) / h
     #z - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=3, Nx-2, 1
      do k=Nz, Nz, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D
     #-1 * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2
     # * np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, 
     #k)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.41666666666666
     #67D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 
     #* n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / h
     #z - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=2, 2, 1
      do k=Nz, Nz, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D
     #-1 * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2
     # * np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, 
     #k)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.41666666666666
     #67D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 
     #* n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / h
     #z - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=1, 1, 1
      do k=Nz, Nz, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D
     #-1 * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2
     # * np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, 
     #k)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.41666666666666
     #67D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 
     #* n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / h
     #z - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      do i=Nx-1, Nx-1, 1
      do k=Nz, Nz, 1
      qb = (-0.1D1 * n_f(i, k) + np1_f(i, k)) / ht - 0.4166666666666667D
     #-1 * (-0.3D1 * np1_f(i, k - 4) + 0.16D2 * np1_f(i, k - 3) - 0.36D2
     # * np1_f(i, k - 2) + 0.48D2 * np1_f(i, k - 1) - 0.25D2 * np1_f(i, 
     #k)) / hz - 0.5000000000000000D0 * np1_f_t(i, k) - 0.41666666666666
     #67D-1 * (-0.3D1 * n_f(i, k - 4) + 0.16D2 * n_f(i, k - 3) - 0.36D2 
     #* n_f(i, k - 2) + 0.48D2 * n_f(i, k - 1) - 0.25D2 * n_f(i, k)) / h
     #z - 0.5000000000000000D0 * n_f_t(i, k)
      res = res + qb**2
      end do
      end do
      res = sqrt(res/(1*Nx*Nz))
      END
