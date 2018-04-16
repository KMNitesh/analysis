c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ect1  --  charge tranfer energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ect1" calculates the charge tranfer energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine ect1
      use sizes
      use atoms
      use ct
      use mpole
      use bound
      use deriv
      use energi
      use group
      use usage
      use virial
      use inform
      use mutant
      implicit none
      integer i,ia,ib
      real*8 e,ecto
      real*8 dt,dt2,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 viro(3,3)
      real*8, allocatable :: decto(:,:)
      logical proceed
      real*8 miu, eta, qct, r1 , r2
      real*8 qct2, qctdr
      logical muti,mutk
c
c
c     zero out the charge tranfer energy and first derivatives
c
      ect = 0.0d0
      do i = 1, n
         dect(1,i) = 0.0d0
         dect(2,i) = 0.0d0
         dect(3,i) = 0.0d0
      end do
c
c     calculate the charge tranfer energy and first derivatives
c
      if (debug) then
         write(*,*) "Charge Tranfer Terms :",nct
      end if

      do i = 1, nct
         ia = ia_ct(i)
         ib = ib_ct(i)
         muti = mut(ia)
         mutk = mut(ib)

c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (proceed)  proceed = (use(ia) .or. use(ib))

         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            call image (xab,yab,zab)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            if (rab < 8) then
                miu = ct_miu(i)
                eta = ct_eta(i)
                qct = ct_qct(i)
                r1  = ct_r1(i)
                r2  = ct_r2(i)
c                if (rab .lt. 2.0d0) then
c                    rab = 2.0d0
c                end if
c                qct2 = qct * exp(-rab * r2)
c                qctdr = - qct2 * r2
                if (rab .gt. r2) cycle
                if (rab .lt. r1) then
                    qct2 = qct
                    qctdr = 0.0d0
                else
                    qct2 = 0.5d0 * qct * (1.0d0 +
     &                 cos(3.1415926d0*(rab-r1)/(r2-r1)))
                    qctdr = -0.5d0 * qct *
     &        sin(3.1415926d0*(rab-r1)/(r2-r1))* 3.1415926d0 /(r2-r1)
                end if

                if ((muti .and. .not.mutk) .or.
     &                (mutk .and. .not.muti)) then
                    e = - miu * abs(qct2) + 0.5 * eta * qct2 * qct2
                    e = e * ctlambda
                    deddt = -miu * qctdr + eta * qct2 * qctdr
                    deddt = deddt * ctlambda
                else
                    e = - miu * abs(qct2) + 0.5 * eta * qct2 * qct2
c                rpole(1,ia) = rpole(1,ia) - qct2
c                rpole(1,ib) = rpole(1,ib) + qct2
                    deddt = -miu * qctdr + eta * qct2 * qctdr
                end if 

                de = deddt / rab
c                if (rab .lt. 2.0d0) then
c                    de = 0.0d0
c                end if
                dedx = de * xab
                dedy = de * yab
                dedz = de * zab
                ect = ect + e
                dect(1,ia) = dect(1,ia) + dedx
                dect(2,ia) = dect(2,ia) + dedy
                dect(3,ia) = dect(3,ia) + dedz
                dect(1,ib) = dect(1,ib) - dedx
                dect(2,ib) = dect(2,ib) - dedy
                dect(3,ib) = dect(3,ib) - dedz
c
c     increment the internal virial tensor components
c
                vxx = xab * dedx
                vyx = yab * dedx
                vzx = zab * dedx
                vyy = yab * dedy
                vzy = zab * dedy
                vzz = zab * dedz
                vir(1,1) = vir(1,1) + vxx
                vir(2,1) = vir(2,1) + vyx
                vir(3,1) = vir(3,1) + vzx
                vir(1,2) = vir(1,2) + vyx
                vir(2,2) = vir(2,2) + vyy
                vir(3,2) = vir(3,2) + vzy
                vir(1,3) = vir(1,3) + vzx
                vir(2,3) = vir(2,3) + vzy
                vir(3,3) = vir(3,3) + vzz
            end if
         end if
      end do

      return
      end
