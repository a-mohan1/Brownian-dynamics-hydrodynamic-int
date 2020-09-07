c     This subroutine uses the semi-implicit predictor-corrector method
c     of Somasi et al., Hsieh et al. and Schroeder et al.
c     to propagate the simulation one step forward in time
c     FREE-DRAINING

c     Inputs:
c     xc, yc, zc = bead coords at current time step istep
c     xn, yn, zn = bead coords at next time step istep+1
c     nbeads = number of beads
c     istep = current time step
c     sprtype = spring force law (no springs = 0, Hookean = 1, Fraenkel = 2, 
c                                 FENE = 3, WLS = 4)
c     Hspr = spring constant
c     whicRNG = which random number generator to use (Knuth = 1)
c     ranflag = 0 at first call to RNG, 1 otherwise; supplied by calling routine
c     nseed = seed for RNG (negative integer for whicRNG = 1)
c     delt = time step
c     Nks = number of Kuhn lengths per spring
c     v = EV parameter, used only if isgev is .TRUE.
c     kappa = 3-by-3 transposed vel gradient tensor (non-dimensionalized with the
c     convective time scale) times the Peclet number, used only if isflow is .TRUE.
c     tether = .TRUE. if there is a wall and the 1st bead is tethered to it (param)
c     overext is set to .TRUE. if subroutine calcFspr encounters an overstretched
c     FENE or WLC spring
c     isgev = .TRUE. if a Gaussian EV force must be included (param)
c     isflow = .TRUE. is there is an imposed (linear) fluid flow (param)
c     istab indicates whether or not the look-up table for a WLS or FENE spring 
c     has been created and saved. It is initialized to .FALSE. in the program
c     and is changed to .TRUE. after the table is created during the first call 
c     to subroutine Qfromrhs; reqd only for FENE and WLC springs

      subroutine pc1sfd(xc, yc, zc, xn, yn, zn, nbeads, istep,
     &                sprtype, Hspr, whicRNG, ranflag, nseed, delt, Nks,
     &                v, kappa, tether, overext, isgev, isflow, istab)

      implicit none
      integer nbeads, istep, sprtype, whicRNG, ranflag, nseed
      integer i, icoord, jcoord, ibead, icount, jcount
      real*8 xc(nbeads), yc(nbeads), zc(nbeads)
      real*8 xn(nbeads), yn(nbeads), zn(nbeads)
      real*8 Hspr, delt, Nks, v, kappa(3, 3), rannum
      real*8 Qc(nbeads-1, 3), Qprime(3), Qc1(nbeads-1, 3)
      real*8 Qf(nbeads-1, 3), const
      real*8 Qvect(3), fspt(3), fspful(3*(nbeads-1))
      real*8 fevful(3*nbeads), evforc(3), n(3*nbeads), tevb(nbeads-1, 3)
      real*8 rhsc1(3), rhsc2(3)
      real*8 res, tol
      real*8 rold(3*nbeads), rnew(3*nbeads)
      logical tether, overext, isgev, isflow, istab
      
c     i = spring index, runs from 1..nbeads-1
c     ibead = bead index, run from 1..nbeads
c     icoord, jcoord = coord indices, run from 1..3
c     icount, jcount = index from 1..3*nbeads
c     (or 1..3*(nbeads-1), if used over springs)

c     Function rannum returns a uniform deviate in [0, 1]

c     Qc = current vector of spring vectors, with Qc(i, 1) being the
c     x-comp, Qc(i, 2) the y-comp and Qc(i, 3) the z-comp (time istep)
c     Qprime = spring vector of spring i after the predictor step
c     Qprime is repeatedly overwritten while iterating over springs i
c     Qc1 = spring vectors after the 1st corrector step 
c     Qf = final spring vectors, after the 2nd corrector step
c     Qvect = spring vector of a specified spring
c     fspt = force vector from a specified spring
c     (Qvect and fspt are temporary vars, repeatedly overwritten)
c     fspful = force vector from all springs initially at current time step
c     (components of any one spring are in contiguous locations)
c     fspful will be overwritten after 1st and 2nd corrector steps
c     as and when Qc1 and Qf become available

c     const = sqrt(24 dt)
      const = sqrt(24.0*delt)
c     res = residual at the end of 2nd corrector step
c     tol = tolerance (1e-6)
      tol = 1.d-6

c     fevful = EV force on all beads, computed once at the start of istep
c     Forces on each bead are in contiguous locations; one bead follows another,
c     so fevful(3*(j-1)+1,2,3) is the EV force vector on bead j
c     evforc = EV force vector on a specified bead
c     (temporary var, repeatedly overwritten)
c     n = 3*nbeads-by-1 vector of uniform rvs in [-0.5, 0.5]
c     computed at the start of istep
c     tevb = stores the EV and Brownian contributions, to be
c     transferred to the appropriate vars rhsc1 and rhsc2

c     rhsc1, rhsc2 = right hand side vectors of the eqs 
c     obtained in the 1st and 2nd corrector steps
c     These vectors are overwritten in each iteration over springs i

c     rold = array of coords of all beads at time istep
c     rnew = array of coords of all beads at the end of predictor step
c     (coords of any bead are in contiguous locations)

c     Assemble the spring vectors Q_i = r_{i+1}-r_{i}
c     from the position vectors of beads i+1 (r_{i+1}) and i (r_i)
c     at the current time point, istep

      do 10 i = 1, nbeads-1
        Qc(i, 1) = xc(i+1) - xc(i)        
        Qc(i, 2) = yc(i+1) - yc(i)
        Qc(i, 3) = zc(i+1) - zc(i)
 10   continue

c     Assemble the 3*nbeads-by-1 EV force vector fevful, if necessary,
c     at time istep
c     or set it to zero if there is no Gaussian EV force 
      if(isgev) then
       do 20 ibead = 1, nbeads
        call Fev1b(xc, yc, zc, nbeads, ibead, v, Nks, evforc)  
        do 30 icoord = 1, 3
         icount = 3*(ibead-1)+icoord
         fevful(icount) = evforc(icoord)
 30     continue
 20    continue
      else
       do 40 icount = 1, 3*nbeads
          fevful(icount) = 0.d0
 40    continue
      endif

c     Now calculate the force from all springs at istep
c     This will be reqd in the predictor and 1st corrector step
c     Always check for overextended springs

c     overext is initialized to .FALSE. in calcFspr 
c     and need not be initialized here
c     In case any over-extended springs are encountered by subroutine
c     calcFspr, the current time step, istep, must be recomputed, 
c     before the simulation can be moved forward in time to istep+1,
c     so simply return control to the calling routine whenever overext
c     takes the value .TRUE. after a call to calcFspr

c     Qvect = spring vector of spring i 
c     fspful = force vector from all springs at current time step
c     (components for any one spring are in contiguous locations)
c     fspt = force vector from spring i

      do 45 i = 1, nbeads - 1
        do 50 icoord = 1, 3
            Qvect(icoord) = Qc(i, icoord)
 50     continue
        call calcFspr(sprtype, Qvect, Hspr, fspt, overext)    
        if(overext) goto 5000
        do 55 icoord = 1, 3
            icount = 3*(i-1)+icoord
            fspful(icount) = fspt(icoord)
 55     continue
 45   continue

c     Now generate a vector n of uniform deviates in [-0.5, 0.5]
c     The same 3nbeads-by-1 vector n must be used
c     throughout this time step
c     Function rannum returns a uniform deviate in [0, 1]
      do 60 icount = 1, 3*nbeads
          n(icount) = rannum(whicRNG, ranflag, nseed) - 0.5
 60   continue

c     The first step is the predictor step - which is the explicit Euler step
c     Bead locations based on the explicit Euler method will be computed 
c     and saved to xn, yn and zn
c     The implicit method is meant to ensure that springs are not overstretched 
c     and will not affect bead 1
c     So the bead 1 location at time istep+1 is obtained from this predictor step

c     BEGINNING OF PREDICTOR STEP
     
c     Assemble bead coords at time istep into 3*nbeads vector rold
      call xyz2r(xc, yc, zc, nbeads, rold)

c     calculate bead coords at time istep+1 (saved in rnew)
c     from a Langevin equation for each coordinate of each bead
c     by the explicit Euler method
  
c     Iterate over beads and coordinates

      do 70 ibead = 1, nbeads
       do 80 icoord = 1, 3
        icount = 3*(ibead-1)+icoord
c       Set rnew to rold, which is the first term in the eq for rnew
        rnew(icount) = rold(icount)

c       Now, if necessary, add the contribution from fluid flow
        if(isflow) then 
         do 85 jcoord = 1, 3
          jcount = 3*(ibead-1)+jcoord
          rnew(icount) = rnew(icount) + 
     &    kappa(icoord, jcoord)*rold(jcount)*delt
 85      continue
        endif

c       Add the spring forces
c       Must treat ibead=1 (first bead) and ibead=nbeads (last bead) separately
c       For ibead=1, spring force is only from spring 1     
c       (first 3 elements of fspful)
c       For ibead=nbeads, spring force is only from spring nbeads-1
c       (last 3 elements of fspful) and takes a negative sign
c       For bead ibead, take the diff between forces from springs ibead and ibead-1

        if(ibead.eq.1) then
            rnew(icount) = rnew(icount) + 
     &      fspful(icount)*delt
        elseif(ibead.eq.nbeads) then                
            rnew(icount) = rnew(icount) + 
     &      (-fspful(icount-3))*delt
        else
            rnew(icount) = rnew(icount) + 
     &      (fspful(icount)-fspful(icount-3))*delt
        endif 

c       Now add the EV force terms, if necessary
        if(isgev) then
           rnew(icount) = rnew(icount) + fevful(icount)*delt
        endif

c       Now add the Brownian terms
        rnew(icount) = rnew(icount) + const*n(icount)

c       This completes the predictor step for bead ibead, coord icoord
c       Finish iterations    
 80    continue
 70   continue

c     If bead 1 is tethered, it must not be moved
      if(tether) then 
c     revert to old bead 1 coordinates
         rnew(1) = rold(1)
         rnew(2) = rold(2)
         rnew(3) = rold(3)
      endif

c     calculate xn, yn, zn using rnew
      call r2xyz(rnew, nbeads, xn, yn, zn) 

c     END OF PREDICTOR STEP
c     We now know x, y, z at the end of an explicit Euler step
      
c     Will compute the EV and Brownian terms for the spring eqs and store them, 
c     since the same terms are used in both corrector steps, throughout istep

      do 110 i = 1, nbeads-1
        do 120 icoord = 1, 3  
c        Add the EV and Brownian terms
         tevb(i, icoord) = delt*
     &   (fevful(3*i+icoord)-fevful(3*(i-1)+icoord))
     &   + const*(n(3*i+icoord)-n(3*(i-1)+icoord))  
 120    continue
 110  continue
         
c     BEGINNING OF FIRST CORRECTOR STEP FOR ALL SPRINGS

c     Iterate over springs
      do 150 i = 1, nbeads-1  

c      First, calculate Qprime for i from the x, y, z values after predictor step
       Qprime(1) = xn(i+1) - xn(i)
       Qprime(2) = yn(i+1) - yn(i)
       Qprime(3) = zn(i+1) - zn(i)
 
c      BEGINNING OF FIRST CORRECTOR STEP FOR SPRING i
c      start iterations over coordinates
       do 200 icoord = 1, 3

c       This step involves the solution of a (possibly cubic) equation for Qc1
c       Begin by calculating the right hand side vector rhsc1
c       of the eq for the current spring i in the current coord direction icoord

c       Start with the terms already known at this point
        rhsc1(icoord) = Qc(i, icoord) + tevb(i, icoord)

c       Now if needed, calculate the flow term
        if(isflow) then
          do 220 jcoord = 1, 3
           rhsc1(icoord) = rhsc1(icoord) + 
     &           0.5*delt*kappa(icoord, jcoord)* 
     &           (Qprime(jcoord)+Qc(i, jcoord))
 220      continue
        endif

c       Now add the spring force terms
c       For springs<i, Qc1 vector is known, and must be used to update fspful
c       For springs>=i, must take the spring forces from fspful at time istep

c       As and when they become available, the fspful vector will be updated
c       with the forces calculated from Qc1

c       Can now use the most up-to-date fspful vector to calculate 
c       the spring force contributions as was done in the predictor step
c       fspful is calculated using Qc1 for springs < i, and Qc for i and beyond

        if(i.eq.1) then
           rhsc1(icoord) = rhsc1(icoord) + fspful(3*i+icoord)*delt
        elseif(i.eq.nbeads-1) then                
           rhsc1(icoord) = rhsc1(icoord) + fspful(3*(i-2)+icoord)*delt
        else
           rhsc1(icoord) = rhsc1(icoord) + 
     &     (fspful(3*i+icoord)+fspful(3*(i-2)+icoord))*delt
        endif 

c       Finish iterations over coordinates
 200   continue
c      The rhsc1 vector is now available for spring i

c      Now solve the relevant equation to get Qc1 for spring i
c      Qc1 vector of spring i is temporarily stored in Qvect
       call Qfromrhs(rhsc1, sprtype, Hspr, delt, Qvect, istab)
c      Qc1 vector for spring i is now available

c      Must now update fspful
c      since Qc1 for spring i is available

       call calcFspr(sprtype, Qvect, Hspr, fspt, overext)
       if(overext) goto 5000
       do 350 jcoord = 1, 3
            jcount = 3*(i-1)+jcoord
            fspful(jcount) = fspt(jcoord)       
 350   continue

c      Now copy Qvect into Qc1
       do 380 jcoord = 1, 3
          Qc1(i, jcoord) = Qvect(jcoord)
 380   continue
  
c      END OF FIRST CORRECTOR STEP FOR SPRING i
c      Qc1 for spring i is now available
c      and fspful has been updated using Qc1 for spring i

c      Finish iterations over springs i
 150  continue 
c     END OF FIRST CORRECTOR STEP FOR ALL SPRINGS

c     We now have Qc1 for all springs 
c     and fspful calculated from Qc1 for all springs


c     THIS IS WHERE WE LAND IF WE FAIL TO REDUCE THE RESIDUAL BELOW THE TOLERANCE
 390  continue

c     BEGINNING OF SECOND CORRECTOR STEP FOR ALL SPRINGS

      do 400 i = 1, nbeads - 1

c      BEGINNING OF SECOND CORRECTOR STEP FOR SPRING i

c      begin iterations over coordinates

       do 410 icoord = 1, 3

c       This step involves the calculation of Qf
c       Begin by calculating the right hand side vector rhsc2
c       of the eq for the current spring i in the current coord direction icoord

c       Start with the terms already known at this point
        rhsc2(icoord) = Qc(i, icoord) + tevb(i, icoord)
         
c       Now if needed, calculate the flow term
        if(isflow) then
          do 420 jcoord = 1, 3
           rhsc2(icoord) = rhsc2(icoord) + 
     &           0.5*delt*kappa(icoord, jcoord)* 
     &           (Qc1(i, jcoord)+Qc(i, jcoord))
 420      continue
        endif

c       Now add the spring force terms
c       For springs<i, Qf vector is known, and must be used to update fspful
c       For springs>=i, must take the spring forces from fspful computed using Qc1

c       As and when they become available, the fspful vector will be updated
c       with the forces calculated from Qf

c       Can now use the most up-to-date fspful vector to calculate 
c       the spring force contributions as was done in the 1st corrector step
c       fspful is calculated using Qf for springs < i, and Qc1 for i and beyond

        if(i.eq.1) then
           rhsc2(icoord) = rhsc2(icoord) + fspful(3*i+icoord)*delt
        elseif(i.eq.nbeads-1) then                
           rhsc2(icoord) = rhsc2(icoord) + fspful(3*(i-2)+icoord)*delt
        else
           rhsc2(icoord) = rhsc2(icoord) + 
     &     (fspful(3*i+icoord)+fspful(3*(i-2)+icoord))*delt
        endif 
 
c       Finish iterations over coordinates
 410   continue
c      The rhsc2 vector is now available for spring i

c      Now solve the relevant equation to get Qf for spring i
c      Qf for spring i is temporarily stored in Qvect
       call Qfromrhs(rhsc2, sprtype, Hspr, delt, Qvect, istab)
c      Qf vector for spring i is now available

c      Must now update fspful
c      since Qf for spring i is available

       call calcFspr(sprtype, Qvect, Hspr, fspt, overext)
       if(overext) goto 5000
       do 500 jcoord = 1, 3
            jcount = 3*(i-1)+jcoord
            fspful(jcount) = fspt(jcoord)       
 500   continue
c  
c      Now transfer Qvect to Qf for i
       do 520 jcoord = 1, 3
          Qf(i, jcoord) = Qvect(jcoord)
 520   continue

c      END OF SECOND CORRECTOR STEP FOR SPRING i
c      Qf for spring i is now available
c      and fspful has been updated
c      Finish iterations over springs 

 400  continue
c     END OF SECOND CORRECTOR STEP FOR ALL SPRINGS

c     We now know Qf for all springs 
c     and fspful calculated from Qf for all springs

c     Now compute the residual and check if it is acceptable
c     If not acceptable, then transfer Qf to Qc1 and goto 390
c     Recall that now fspful contains values computed from Qf 
c     (i.e., the new Qc1)

      res = 0
      do 550 i = 1, nbeads - 1
         res = res + (Qf(i, 1) - Qc1(i, 1))**2 + 
     &               (Qf(i, 2) - Qc1(i, 2))**2 + 
     &               (Qf(i, 3) - Qc1(i, 3))**2
 550  continue     
      res = sqrt(res)
      
      if(res.gt.tol) then
        do 600 i = 1, nbeads - 1
         do 620 icoord = 1, 3 
          Qc1(i, icoord) = Qf(i, icoord)
 620     continue
 600    continue
        goto 390
      endif 

c     We end up here if residual <= tolerance
c     Must now compute xn, yn and zn at istep + 1

c     Bead 1 coordinates were already known after the predictor step
c     (tethered or otherwise)
    
c     Can now compute coords of beads 2 to nbeads by adding spring vector
c     of connecting spring to preceding bead coords 

      do 640 ibead = 2, nbeads
       xn(ibead) = xn(ibead-1) + Qf(ibead-1, 1) 
       yn(ibead) = yn(ibead-1) + Qf(ibead-1, 2) 
       zn(ibead) = zn(ibead-1) + Qf(ibead-1, 3)
 640  continue 

   
5000  continue

      return
      end

             
