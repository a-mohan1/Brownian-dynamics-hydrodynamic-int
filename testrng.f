c     This program tests function rannum

      program testrng

      implicit none
      integer whicRNG, ranflag, nseed, i, num
      real*8 rannum, rvs(10000), av, var

c     initialize ranflag to 0 for first call to RNG
      ranflag = 0

c     number of uniform rvs generated (max. 10 000)
      num = 10000

c     initialize mean and variance
c     uniform distr in [0, 1] has mean 0.5, variance 1/12 = 0.083
      av = 0.d0
      var = 0.d0 

c     whicRNG = 0 if ranuls, 1 if ran3
c     if ranuls, nseed must be an integer between 0 and 2E9
c     if ran3, nseed must be a negative integer
      whicRNG = 1

c     nseed must now be a negative integer 
      nseed = -2006
 
      do 10 i = 1, num
         rvs(i) = rannum(whicRNG, ranflag, nseed)
         write(*, *) rvs(i)
         av = av + rvs(i)
         var = var + rvs(i)**2
 10   continue

      av = av/num
      var = var/num - av**2  
      write(*, *) av, var

      stop
      end




