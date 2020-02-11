c
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kumbrella  --  US parameter assignment      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kumbella" assigns umbella sampling parameters and options
c
c
      subroutine kumbrella
      use sizes
      use atoms
      use bound
      use boxes
      use chunks
      use umbrella
      use inform
      use iounit
      use keys
      use limits
      use openmp
      implicit none
      integer next,i,j,k
      integer list(20)
      character*20 keyword
      character*120 record
      character*120 string
c
c     search keywords for Ewald summation commands
c
      pull_switch = 0
      pull_group1_num = 0
      pull_group2_num = 0
      current_step = 0
      do i = 1, nkey
         record = keyline(i)
         next = 1
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:120)
         if (keyword(1:5) .eq. 'PULL ') then
             pull_switch = 1
         else if (keyword(1:11) .eq. 'PULL-CONST ') then
            read (string,*,err=20) pull_const
         else if (keyword(1:10) .eq. 'PULL-DIST ') then
            read (string,*,err=20) pull_dist
         else if (keyword(1:15) .eq. 'PULL-DIST-RATE ') then
            read (string,*,err=20) pull_dist_rate
         else if (keyword(1:12) .eq. 'PULL-GROUP1 ') then
            string = record(next:120)
            do k = 1, 20
                list(k) = 0
            end do
            read (string,*,err=10,end=10) (list(k),k=1,20)
   10       continue
            k = 1
            do while (list(k) .ne. 0)
                if (list(k) .gt. 0) then
                    j = list(k)
                    pull_group1_num = pull_group1_num + 1
                    pull_group1(pull_group1_num) = j
                    k = k + 1
                else
                    do j = abs(list(k)), abs(list(k+1))
                        pull_group1_num = pull_group1_num + 1
                        pull_group1(pull_group1_num) = j
                    end do
                    k = k + 2
                end if
            end do
         else if (keyword(1:12) .eq. 'PULL-GROUP2 ') then
            string = record(next:120)
            do k = 1, 20
                list(k) = 0
            end do
            read (string,*,err=30,end=30) (list(k),k=1,20)
   30       continue
            k = 1
            do while (list(k) .ne. 0)
                if (list(k) .gt. 0) then
                    j = list(k)
                    pull_group2_num = pull_group2_num + 1
                    pull_group2(pull_group2_num) = j
                    k = k + 1
                else
                    do j = abs(list(k)), abs(list(k+1))
                        pull_group2_num = pull_group2_num + 1
                        pull_group2(pull_group2_num) = j
                    end do
                    k = k + 2
                end if
            end do
         end if
   20    continue
      end do
      if (verbose) then
         write (iout,*) pull_const,pull_dist,pull_dist_rate
      end if
      return
      end
