      subroutine u_pp(a,x,alpha,n_pi,n_pp,np1_pi,np1_pp,Nx,ht,hx,epsdis,
     &myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 epsdis
      real*8 myzero
      real*8 a(Nx)
      real*8 x(Nx)
      real*8 alpha(Nx)
      real*8 n_pi(Nx)
      real*8 n_pp(Nx)
      real*8 np1_pi(Nx)
      real*8 np1_pp(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      include 'tvd.inc'
      do i=1, 1, 1
      qb = -0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      do i=2, 2, 1
      qb = np1_pp(i) - 0.1D1 * ht * (-0.1D1 * (n_pp(i) - 0.1D1 * np1_pp(
     #i)) / ht + 0.2500000000000000D0 * (alpha(i - 1) - 0.1D1 * alpha(i 
     #+ 1)) / hx / a(i) * np1_pi(i) - 0.2500000000000000D0 * alpha(i) / 
     #a(i) ** 2 * np1_pi(i) * (a(i - 1) - 0.1D1 * a(i + 1)) / hx - 0.250
     #0000000000000D0 * alpha(i) / a(i) * (-0.1D1 * np1_pi(i - 1) + np1_
     #pi(i + 1)) / hx + 0.2500000000000000D0 * (alpha(i - 1) - 0.1D1 * a
     #lpha(i + 1)) / hx / a(i) * n_pi(i) - 0.2500000000000000D0 * alpha(
     #i) / a(i) ** 2 * n_pi(i) * (a(i - 1) - 0.1D1 * a(i + 1)) / hx - 0.
     #2500000000000000D0 * alpha(i) / a(i) * (-0.1D1 * n_pi(i - 1) + n_p
     #i(i + 1)) / hx)
      res(i)=qb
      end do
      do i=3, Nx-2, 1
      qb = np1_pp(i) - 0.16D2 / (0.16D2 + 0.3D1 * epsdis) * ht * (-0.1D1
     # * (n_pp(i) - 0.1D1 * np1_pp(i)) / ht + 0.2500000000000000D0 * (al
     #pha(i - 1) - 0.1D1 * alpha(i + 1)) / hx / a(i) * np1_pi(i) - 0.250
     #0000000000000D0 * alpha(i) / a(i) ** 2 * np1_pi(i) * (a(i - 1) - 0
     #.1D1 * a(i + 1)) / hx - 0.2500000000000000D0 * alpha(i) / a(i) * (
     #-0.1D1 * np1_pi(i - 1) + np1_pi(i + 1)) / hx + 0.3125000000000000D
     #-1 * epsdis / ht * (0.6D1 * np1_pp(i) + np1_pp(i + 2) + np1_pp(i -
     # 2) - 0.4D1 * np1_pp(i + 1) - 0.4D1 * np1_pp(i - 1)) + 0.250000000
     #0000000D0 * (alpha(i - 1) - 0.1D1 * alpha(i + 1)) / hx / a(i) * n_
     #pi(i) - 0.2500000000000000D0 * alpha(i) / a(i) ** 2 * n_pi(i) * (a
     #(i - 1) - 0.1D1 * a(i + 1)) / hx - 0.2500000000000000D0 * alpha(i)
     # / a(i) * (-0.1D1 * n_pi(i - 1) + n_pi(i + 1)) / hx + 0.3125000000
     #000000D-1 * epsdis / ht * (0.6D1 * n_pp(i) + n_pp(i + 2) + n_pp(i 
     #- 2) - 0.4D1 * n_pp(i + 1) - 0.4D1 * n_pp(i - 1)))
      res(i)=qb
      end do
      do i=Nx-1, Nx-1, 1
      qb = np1_pp(i) - 0.1D1 * ht * (-0.1D1 * (n_pp(i) - 0.1D1 * np1_pp(
     #i)) / ht + 0.2500000000000000D0 * (alpha(i - 1) - 0.1D1 * alpha(i 
     #+ 1)) / hx / a(i) * np1_pi(i) - 0.2500000000000000D0 * alpha(i) / 
     #a(i) ** 2 * np1_pi(i) * (a(i - 1) - 0.1D1 * a(i + 1)) / hx - 0.250
     #0000000000000D0 * alpha(i) / a(i) * (-0.1D1 * np1_pi(i - 1) + np1_
     #pi(i + 1)) / hx + 0.2500000000000000D0 * (alpha(i - 1) - 0.1D1 * a
     #lpha(i + 1)) / hx / a(i) * n_pi(i) - 0.2500000000000000D0 * alpha(
     #i) / a(i) ** 2 * n_pi(i) * (a(i - 1) - 0.1D1 * a(i + 1)) / hx - 0.
     #2500000000000000D0 * alpha(i) / a(i) * (-0.1D1 * n_pi(i - 1) + n_p
     #i(i + 1)) / hx)
      res(i)=qb
      end do
      do i=Nx, Nx, 1
      qb = np1_pp(i) - 0.4D1 / (0.4D1 * hx * x(i) + 0.3D1 * ht * x(i) + 
     #0.2D1 * ht * hx) * ht * hx * x(i) * (-0.1D1 * (n_pp(i) - 0.1D1 * n
     #p1_pp(i)) / ht - 0.2500000000000000D0 * (-0.1D1 * np1_pp(i - 2) + 
     #0.4D1 * np1_pp(i - 1) - 0.3D1 * np1_pp(i)) / hx + 0.50000000000000
     #00D0 * np1_pp(i) / x(i) - 0.2500000000000000D0 * (-0.1D1 * n_pp(i 
     #- 2) + 0.4D1 * n_pp(i - 1) - 0.3D1 * n_pp(i)) / hx + 0.50000000000
     #00000D0 * n_pp(i) / x(i))
      res(i)=qb
      end do
      END
