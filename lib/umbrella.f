c
c     ################################################################
c     ##                                                            ##
c     ##  module umbrella  --  Umbrella Sampling options            ##
c     ##                                                            ##
c     ################################################################
c
c
c     pull_const   force constant of harmonic potential(kcal/mol/ang^2
c     pull_dist    distance to maintain in US windows
c     pull_dist_rate  pull_dist linear increase rate (ang/step)
c
c
      module umbrella
      use sizes
      implicit none
      real*8 pull_const
      real*8 pull_dist
      real*8 pull_dist_rate
      integer pull_group1(maxatm)
      integer pull_group1_num
      integer pull_group2(maxatm)
      integer pull_group2_num

      integer current_step
      integer pull_switch
      save
      end
