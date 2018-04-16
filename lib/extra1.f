c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra1  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra1" calculates any additional user defined potential
c     energy contribution and its first derivatives
c
c
      subroutine extra1
      use sizes
      use atoms
      use atomid
      use deriv
      use energi
      use umbrella
      use bound
      use iounit
      use virial
      implicit none
      integer i,j,k
      real*8 weigha, weighb, weigh
      real*8 xcm, ycm,zcm
      real*8 xr,yr,zr
      real*8 r,dt,dt2
      real*8 e,de,dedx,dedy,dedz
      real*8 vxx,vyx,vzx,vyy,vzy,vzz
      real*8 ratio
c
c
c     zero out the extra energy term and first derivatives
c
      ex = 0.0d0
      do i = 1, n
         dex(1,i) = 0.0d0
         dex(2,i) = 0.0d0
         dex(3,i) = 0.0d0
      end do
c
c     add any user-defined extra potentials and derivatives;
c     also increment intermolecular energy and virial as needed
c
c     e = ......
c     ex = ex + e
c     do i = 1, n
c        dex(1,i) = ......
c        dex(2,i) = ......
c        dex(3,i) = ......
c     end do
c
c     ADD  Umbrella Sampling Code
      if ( pull_switch .eq. 0) then
        return
      end if
c     Calculate the center of mass of group1

      weigha = 0.0
      xcm = 0.0
      ycm = 0.0
      zcm = 0.0
      do j = 1, pull_group1_num
        k = pull_group1(j)
        weigh = mass(k)
        weigha = weigha + weigh
        xcm = xcm + x(k) * weigh
        ycm = ycm + y(k) * weigh
        zcm = zcm + z(k) * weigh

      end do
      xr = xcm / weigha
      yr = ycm / weigha
      zr = zcm / weigha

      weighb = 0.0
      xcm = 0.0
      ycm = 0.0
      zcm = 0.0
      do j = 1, pull_group2_num
        k = pull_group2(j)
        weigh = mass(k)
        weighb = weighb + weigh
        xcm = xcm + x(k) * weigh
        ycm = ycm + y(k) * weigh
        zcm = zcm + z(k) * weigh
      end do
      xr = xr - xcm/weighb
      yr = yr - ycm/weighb
      zr = zr - zcm/weighb

      if (use_bounds)  call image (xr,yr,zr)

      r = sqrt(xr*xr+yr*yr+zr*zr)
      dt = r - (pull_dist + pull_dist_rate *  current_step)
      dt2 = dt * dt
      e = 0.5 * pull_const * dt2

      de = pull_const * dt / r
      dedx = de * xr
      dedy = de * yr
      dedz = de * zr

      ex = ex + e

      do j = 1, pull_group1_num
        k = pull_group1(j)
        ratio = mass(k) / weigha
        dex(1,k) = dex(1,k) + dedx * ratio
        dex(2,k) = dex(2,k) + dedy * ratio
        dex(3,k) = dex(3,k) + dedz * ratio
      end do
      do j = 1, pull_group2_num
        k = pull_group2(j)
        ratio = mass(k) / weighb
        dex(1,k) = dex(1,k) - dedx * ratio
        dex(2,k) = dex(2,k) - dedy * ratio
        dex(3,k) = dex(3,k) - dedz * ratio
      end do
c
c     increment the internal virial tensor components
c
      vxx = xr * dedx
      vyx = yr * dedx
      vzx = zr * dedx
      vyy = yr * dedy
      vzy = zr * dedy
      vzz = zr * dedz
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz

      return
      end
