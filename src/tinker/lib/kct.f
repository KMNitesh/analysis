c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kct  --  charge tranfer term assignment       ##
c     ##                                                           ##
c     ###############################################################
c
c
c
c
      subroutine kct
      use sizes
      use atomid
      use atoms
      use ct
      use couple
      use fields
      use inform
      use iounit
      use kbonds
      use keys
      use potent
      use usage
      implicit none
      integer i,j,t
      integer ia_type,ia_type_t
      integer ib_type,ib_type_t

      nct = 0

c      write(*,*) "nctp:",nctp
      if (nctp.eq.0) then
          use_ct = .false.
          return
      end if

      if (.not. allocated(ia_ct)) allocate (ia_ct(int(n*n/2)))
      if (.not. allocated(ib_ct)) allocate (ib_ct(int(n*n/2)))
      if (.not. allocated(ct_miu))  allocate (ct_miu(int(n*n/2)))
      if (.not. allocated(ct_eta))  allocate (ct_eta(int(n*n/2)))
      if (.not. allocated(ct_qct))  allocate (ct_qct(int(n*n/2)))
      if (.not. allocated(ct_r1))  allocate (ct_r1(int(n*n/2)))
      if (.not. allocated(ct_r2)) allocate (ct_r2(int(n*n/2)))

      do i = 1, n-1
        do j = i+1, n
            ia_type = type(i)
            ib_type = type(j)
            do t = 1, nctp
                ia_type_t = ia_ctp(t)
                ib_type_t = ib_ctp(t)
                if (ia_type_t .eq. ia_type .and.
     &                 ib_type_t .eq. ib_type) then
                    nct = nct +1 
                    ia_ct(nct) = i
                    ib_ct(nct) = j
                    ct_miu(nct) = ct_miup(t)
                    ct_eta(nct) = ct_etap(t)
                    ct_qct(nct) = ct_qctp(t)
                    ct_r1(nct) = ct_r1p(t)
                    ct_r2(nct) = ct_r2p(t)
                else if (ia_type_t .eq. ib_type .and.
     &                 ib_type_t .eq. ia_type) then
                    nct = nct +1 
                    ia_ct(nct) = j
                    ib_ct(nct) = i
                    ct_miu(nct) = ct_miup(t)
                    ct_eta(nct) = ct_etap(t)
                    ct_qct(nct) = ct_qctp(t)
                    ct_r1(nct) = ct_r1p(t)
                    ct_r2(nct) = ct_r2p(t)
                end if
            end do
        end do
      end do
c      write(*,*) "nct:",nct
      if (nct.eq.0) use_ct = .false.
      return
      end


                

