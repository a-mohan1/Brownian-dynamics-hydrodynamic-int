c     This program tests function gauss

      program testgauss

      implicit none
      integer whicRNG, ranflag, jseed, i, num
      real*8 gauss, gdev(10000), av, var

      parameter (whicRNG = 1, jseed = -6)

c     initialize ranflag to 0 for first call to rannum
      ranflag = 0

c     number of rvs generated (max. 10 000)
      num = 10000

c     initialize mean and variance 
c     std normal distr has mean 0, var 1
      av = 0.d0
      var = 0.d0
      
      do 10 i = 1, num
         gdev(i) = gauss(whicRNG, ranflag, jseed)
         write(*, *) gdev(i)
         av = av + gdev(i)
         var = var + gdev(i)*gdev(i)
 10   continue

      av = av/num
      var = var/num - av*av
      write(*, *) av, var

      stop
      end 
