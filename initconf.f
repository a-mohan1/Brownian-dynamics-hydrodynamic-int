c     This program generates initial configurations

c     Variables:
c     xc, yc, zc, xn, yn, zn = bead coords at successive time steps
c     nbeads = number of beads (max. 100)
c     ntsteps = number of time steps to equilibration (max 100 000)
c        since for 100 beads, with H=3Nks/lambda=54, relt=1/(8Hsin^2(pi/2nbeads))
c        =9.38, 5 rel times = 5*9.38/(delt=5e-4), approx 100 000 steps
c     nens = number of ensembles
c     sprtype = spring force law (no springs=0, Hookean=1, Fraenkel=2, 
c                                 FENE=3, WLS=4)
c     lambda = ratio of effective to true persistence length
c     Hspr = spring constant = 3*Nks/lambda
c     whicRNG = which RNG to use to generate uniform deviates (Knuth = 1)
c     ranflag = 0 at first call to RNG, 1 otherwise
c     nseed = seed for RNG (negative integer for whicRNG = 1)
c     delt = time step
c     relt = Rouse relaxation time
c     Nks = number of Kuhn steps per spring
c     v = EV parameter, used if isgev is .TRUE.
c     kappa = 3-by-3 transposed vel gradient tensor (non-dimensionalized by 
c     convective time scale) times Peclet number, not relevant for eqm calc.,
c     used only if there is a linear fluid flow
c     tether = .TRUE. if the first bead is tethered
c     overext is set to .TRUE. if subroutine calcFspr encounters 
c     an overextended FENE spring or WLS
c     isgev is .TRUE. if a Gaussian EV force is included
c     isflow is .TRUE. if there is an imposed linear fluid flow
c     istab indicates whether or not the look-up table for FENE and WLS
c     has been created and saved
c     ibead, istep, ines = counters for beads/ time steps/ ensembles
c     isinp, isconf, isrg indicate whether input and output dat files exist
c     inpstr = temporary character var used to read input string from file
c     Rg = rad. of gyration
c     sampt = time intervals at which sampling is done
c     sstep = time steps after which sampling is done
c     stdev = std dev of spring end-to-end x, y, z coords
c     xinit, yinit, zinit = starting config for all ensembles (gaussian coil)
c     xcm, ycm, zcm = center of mass coords
c     xeq, yeq, zeq = final eqm config
c     icount = counter for time sampling
c     i, j = counters 

c     run for 5 relaxation times
c     Rg sampled at intervals of 0.05 relaxation times
c     xeq, yeq, zeq taken at final time step for last ensemble 
c     (any one realization)

c     all eqs and vars are non-dimensional 

      program initconf

      implicit none
      integer nbeads, ntsteps, nens, sprtype, whicRNG, ranflag, nseed
      integer ibead, istep, iens, sstep, icount, i, j
      character inpstr*10
      real*8 xc(100), yc(100), zc(100)
      real*8 xn(100), yn(100), zn(100)
      real*8 lambda, Hspr, delt, relt, Nks, v, kappa(3, 3)
      real*8 gauss, Rg(100000), sampt, stdev
      real*8 xinit(100), yinit(100), zinit(100)
      real*8 xcm, ycm, zcm, xeq(100), yeq(100), zeq(100)
      logical tether, isgev, isflow, overext, istab, isinp, isconf, isrg

      write(*, *) 'welcome'

c     check if input and output files exist
      inquire(file = 'icdata.dat', exist = isinp)
      if(.NOT.isinp) then
         write(*, *) 'create input file and restart'
         write(*, *) 'bye'
         stop
      endif
      inquire(file = 'iconf.dat', exist = isconf)
      inquire(file = 'rg.dat', exist = isrg)
      if(isconf.OR.isrg) then
         write(*, *) 'delete output file(s) and restart'
         write(*, *) 'bye'
         stop
      endif

c     read input data
      open(unit = 11, file = 'icdata.dat', status = 'OLD')
 10   read(unit = 11, fmt = *, end = 20) inpstr
        backspace(unit=11)
        if(inpstr.eq.'nbeads') then 
          read(11, *) inpstr, nbeads
          write(*, *) 'nbeads = ', nbeads
        elseif(inpstr.eq.'nens') then
          read(11, *) inpstr, nens
          write(*, *) 'ens= ', nens
        elseif(inpstr.eq.'sprtype') then 
          read(11, *) inpstr, sprtype
          write(*, *) 'sprtype= ', sprtype
        elseif(inpstr.eq.'lambda') then
          read(11, *) inpstr, lambda
          write(*, *) 'lambda= ', lambda
        elseif(inpstr.eq.'whicRNG') then
          read(11, *) inpstr, whicRNG
          write(*, *) 'RNG= ', whicRNG
        elseif(inpstr.eq.'nseed') then
          read(11, *) inpstr, nseed
          write(*, *) 'seed= ', nseed 
        elseif(inpstr.eq.'delt') then
          read(11, *) inpstr, delt
          write(*, *) 'time step size= ', delt
        elseif(inpstr.eq.'relt') then
          read(11, *) inpstr, relt
          write(*, *) 'relx time= ', relt
        elseif(inpstr.eq.'Nks') then
          read(11, *) inpstr, Nks
          write(*, *) 'Nks= ', Nks 
        elseif(inpstr.eq.'v') then
          read(11, *) inpstr, v
          write(*, *) 'v= ', v
        elseif(inpstr.eq.'tether') then
          read(11, *) inpstr, tether
          write(*, *) 'tether= ', tether
        elseif(inpstr.eq.'isgev') then
          read(11, *) inpstr, isgev
          write(*, *) 'gEV= ', isgev
        elseif(inpstr.eq.'isflow') then
          read(11, *) inpstr, isflow
          write(*, *) 'linear flow= ', isflow
        else
           write(*, *) 'unrecognized input; bye'
           stop
        endif
        goto 10
 20   continue
      close(unit = 11, status = 'KEEP')

c     define spring constant
      Hspr = 3.d0*Nks/lambda
      write(*, *) 'spring const= ', Hspr

c     run simulation for 5 rel times
      ntsteps = 5.d0*relt/delt
      write(*, *) 'time steps= ', ntsteps

c     sample at intervals of 0.05 rel times
      sampt = 0.05d0*relt
      sstep = sampt/delt
 
c     initialize vars  
      ranflag = 0
      overext = .FALSE.
      istab = .FALSE.
      do 30 icount = 1, ntsteps
         Rg(icount) = 0.d0
 30   continue
      do 40 ibead = 1, nbeads
         xeq(ibead) = 0.d0
         yeq(ibead) = 0.d0
         zeq(ibead) = 0.d0
 40   continue
      
c     define vel grad transpose for linear flow
      do 60 i = 1, 3
       do 50 j = 1, 3
          kappa(i, j) = 0.d0
 50    continue
 60   continue

c     standard dev. of gaussian spring x, y, z end-to-end coords
      stdev = 1.d0/sqrt(3.d0*Nks)
c     generate gaussian coil as starting config for all ensembles
c     function gauss returns a std normal deviate
      xinit(1) = 0.d0
      yinit(1) = 0.d0
      zinit(1) = 0.d0
      do 65 ibead = 2, nbeads
         xinit(ibead) = xinit(ibead-1) +   
     &   gauss(whicRNG, ranflag, nseed)*stdev
         yinit(ibead) = yinit(ibead-1) + 
     &   gauss(whicRNG, ranflag, nseed)*stdev
         zinit(ibead) = zinit(ibead-1) + 
     &   gauss(whicRNG, ranflag, nseed)*stdev
 65   continue

c     run FD simulation   

c     iterate over ensembles
      do 500 iens = 1, nens
         write(*, *) 'ensemble ', iens

c        initialize gaussian coil as starting config for current ensemble
         do 70 ibead = 1, nbeads      
            xc(ibead) = xinit(ibead) 
            yc(ibead) = yinit(ibead) 
            zc(ibead) = zinit(ibead)
 70      continue   

c        initialize counter for current ensemble
         icount = 0

c        iterate over time steps for current ensemble
c        calc bead coords xn, yn, zn at istep from xc, yc, zc at istep-1
         do 400 istep = 1, ntsteps

c           time-stepping
            call pc1sfd(xc, yc, zc, xn, yn, zn, nbeads, istep-1, 
     &      sprtype, Hspr, whicRNG, ranflag, nseed, delt, Nks, v,
     &      kappa, tether, overext, isgev, isflow, istab)

c           check for overextended springs
            if(overext) then
              write(*, *) 'overextended spring encountered; bye'
              stop
            endif

c           if this is a time point at which to sample, calc Rg
            if(mod(istep, sstep).eq.0) then
             icount = icount+1
c            calc center-of-mass
             xcm = xn(1)
             ycm = yn(1) 
             zcm = zn(1)
             do 80 ibead = 2, nbeads
                xcm = xcm + xn(ibead)
                ycm = ycm + yn(ibead)
                zcm = zcm + zn(ibead)
 80          continue
             xcm = xcm/nbeads
             ycm = ycm/nbeads
             zcm = zcm/nbeads
c            calc Rg   
             do 100 ibead = 1, nbeads
               Rg(icount) = Rg(icount) + 
     &         (xn(ibead)-xcm)**2 + (yn(ibead)-ycm)**2 + 
     &         (zn(ibead)-zcm)**2
 100         continue
            endif

c           transfer xn, yn, zn to xc, yc, zc for next time step
            do 150 ibead = 1, nbeads
             xc(ibead) = xn(ibead)
             yc(ibead) = yn(ibead)
             zc(ibead) = zn(ibead)
 150        continue

c           finish iterations over time steps for current ensemble
 400     continue

c        finish iterations over ensembles
 500  continue

c     save config at final time step for last ensemble
c     (we need any one realization)
      do 510 ibead = 1, nbeads
            xeq(ibead) = xn(ibead)
            yeq(ibead) = yn(ibead)
            zeq(ibead) = zn(ibead) 
 510  continue

c     average over ensembles
      do 520 i = 1, icount
         Rg(i) = Rg(i)/(nens*nbeads)
         Rg(i) = sqrt(Rg(i))
 520  continue
        
c     display number of sampled time points
      write(*, *) '#time points at which Rg is sampled= ', icount

c     write to files
      open(unit=12, file='rg.dat', status = 'NEW')
      write(12, 1001) (Rg(i), i = 1, icount)
      open(unit=13, file='iconf.dat', status='NEW')
      do 600 ibead = 1, nbeads
         write(13, 1002) xeq(ibead), yeq(ibead), zeq(ibead) 
 600  continue
      close(unit=12, status='KEEP')
      close(unit=13, status='KEEP')
 1001 format(1e12.6)
 1002 format(3(1X, 1e12.6))
      
      write(*, *) 'op written, bye'

      stop 
      end
