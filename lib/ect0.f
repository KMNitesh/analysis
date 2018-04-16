c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ect0  --  charge tranfer energy                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ect0" calculates the charge tranfer energy
c
c
      subroutine ect0
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
      real*8 e
      real*8 xab,yab,zab,rab
      logical proceed
      real*8 miu, eta, qct, r1 , r2
      real*8 qct2, qctdr
      logical muti,mutk
c
c
c     zero out the charge tranfer energy
c
      ect = 0.0d0
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
                else
                    e = - miu * abs(qct2) + 0.5 * eta * qct2 * qct2
c                rpole(1,ia) = rpole(1,ia) - qct2
c                rpole(1,ib) = rpole(1,ib) + qct2
                end if 
                ect = ect + e
            end if
         end if
      end do

      return
      end
