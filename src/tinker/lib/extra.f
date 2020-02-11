c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine extra  --  user defined extra potentials  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "extra" calculates any additional user defined potential
c     energy contribution
c
c
      subroutine extra
      use sizes
      use atoms
      use atomid
      use energi
      use umbrella
      use bound
      use energi
      implicit none
      integer i,j,k
      real*8 weigha, weighb, weigh
      real*8 xcm, ycm,zcm
      real*8 xr,yr,zr
      real*8 r,dt,dt2
      real*8 e

c
c
c     zero out the energy due to extra potential terms
c
      ex = 0.0d0
c
c     add any user-defined extra potentials below here
c
c     e = ......
c     ex = ex + e
c    ADD  Umbrella Sampling Code
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
      ex = ex + e

      return
      end
