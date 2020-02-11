c
c
c     ################################################################
c     ##                                                            ##
c     ##  module ct  --  charge tranfer term for AMOEBA force field ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxct     maxium number of charge transfer pair in the system
c     nct       total number of charge transfer pair in the system
c     ia_ct     number of the first atom for each pair
c     ib_ct     number of the second atom for each pair
c     ct_miu
c     ct_eta
c     ct_qct
c     ct_r1
c     ct_r2
c
c
      module ct
      implicit none
      integer maxct
      integer nct
      integer, allocatable :: ia_ct(:)
      integer, allocatable :: ib_ct(:)
      real*8, allocatable :: ct_miu(:)
      real*8, allocatable :: ct_eta(:)
      real*8, allocatable :: ct_qct(:)
      real*8, allocatable :: ct_r1(:)
      real*8, allocatable :: ct_r2(:)

      integer maxctp
      parameter (maxctp=1024)
      integer nctp
      integer ia_ctp(maxctp)
      integer ib_ctp(maxctp)
      real*8 ct_miup(maxctp)
      real*8 ct_etap(maxctp)
      real*8 ct_qctp(maxctp)
      real*8 ct_r1p(maxctp)
      real*8 ct_r2p(maxctp)

      save
      end
