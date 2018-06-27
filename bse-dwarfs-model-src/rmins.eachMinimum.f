c Track Rmin and write out 
       subroutine trackRmin(XX,VX,time,newreg)
       INCLUDE "ARCCOM2e2.CH"

       REAL*8 XX(*), VX(*), time
       INTEGER N, iother(3), ll
       LOGICAL newreg

       REAL*8 R2,rminfit,rmincalc,rmincache,rminbest,rminmean,rmb2
       CHARACTER*40 fn

       include "rmins.h"
       logical DOCACHE,transform,coll,implementFlyby

       data iother /3,2,1/

c
c         write (*,*) (XX(i),i=1,3*N)
c         stop

       DOCACHE = .false.
       transform = .false.
       k = 0
       kk = 0
       l = 0
       do i = 1, N-1
          do j = i+1, N
             k = k + 1
             R2 = 0.d0
             do l = 1, 3
                R2 = R2 + (XX(3*i-3+l)-XX(3*j-3+l))**2
             end do
c Test for still approaching minimum
             if (R2.LT.rmin(k)) THEN
                RMIN(k) = R2
                printed(k) = .FALSE.
                docache = .TRUE.
c Test for moving away from a new minimum (so needs "printing")
             else if (.NOT.printed(k)) THEN
c Test for fitting two-body solution to obtain real r_min
                call fitRmin(i,j,k,XX,VX,rminFit,eccfit,tmintau)
c                rminCache = fitRmin(i,j,k,cache(1,1),cache(1,2))
                rminCalc = sqrt(rmin(k))
c                write (22,*) 'Rmins: Fit Calc dt', rminFit, rmincalc, tmintau
                if (rminCalc.lt.rminFit.or.rminCalc.gt.1.d1) then
                   rminBest = rminCalc
                else if (eccfit.lt.1.d-3) then
                   rminBest = rminCalc
                else
                   rminBest = rminFit
c Set a flag so that we know that we need to transform back
                   transform = .true.
c Quick sanity check
                   if (tmintau.gt.3*(time-tcache).or.tmintau.lt.0.d0) then
                      write(*,*) 'WARNING: trackrmin: confused by times'
     : , tmintau, time, tcache, rminCalc, rminFit, eccfit
c                      STOP
                   end if

                end if
                rmb2 = rminBest**2
                rmin(k) = rmb2
                   
                
c Keep global minimum always
                if (rmb2.lt.rminglobal(k)) then
                   rminglobal(k) = rmb2
                   tminglobal(k) = time
                   do ll = 1, 2
                      do l = 1, 3*N
                         cacheglobal(l,ll,k) = cache(l,ll)
                      end do
                   end do
c Store the masses as they can change
                   do l = 1, 3
                      mcache(l) = m(l)
                   end do
                end if

c Implement encounter for minima that are close and non-circular
                if (rmb2.lt.rtid(k).and.eccfit.gt.2.d-3) then
                   WRITE(fn,'("rmin.", i6.6, ".", i1, ".out")')  irun, k
                   OPEN (81,file=fn,access="append")
                   write(81,99001) time,i,j,rminBest,rminCalc,
     :                          (cache(l,1),l=1,3*N),
     :                          (cache(l,2),l=1,3*N),
     :                          (m(l),l=1,N)
99001 FORMAT(1p,e13.5,0p,2(1x,i2),1x,1p,75(e13.5,1x))
                   CLOSE(81)
                   
c Treat tidal capture vs. collision
                   if (rmb2.lt.rcoll(k)) then
                      call implementCollision(k,1)
c Stop this integration now and move on to the next one
                      term = .TRUE.
c Write outcome to statistics file
                      open (80,file="stats.out",access="append")
                      write (80,*) irun,4,i,j,iother(k),
     :                          rmb2,sqrt(rcoll(k))
                      close(80)
                      call flush(80)
                   else
                      coll = implementFlyby(k,time,newreg,xx,vx,tmintau)
                      if (coll) then
                         write (*,*) 'Merger after flyby', irun
c Stop this integration now and move on to the next one
                         term = .TRUE.
c Write outcome to statistics file
                         open (80,file="stats.out",access="append")
                         write (80,*) irun,4,i,j,iother(k),
     :                             rmb2,sqrt(rcoll(k))
                         close(80)
                         call flush(80)
                         call implementCollision(k,2)
                      end if
c                   else if (eccfit.gt.1.d-3) then
c                      call 
c                   else
c                      write (*,*) 'Circular orbit: ignoring minimum', 
c     :                               rminBest, eccfit
                   end if
                end if
                printed(k) = .TRUE.
c If we are finishing the integration here then there is no point completing the loop
                if (term) goto 999

c Moving away from minimum that is already printed; increase "rmin" to catch subsequent minima
             else
c Geometric mean with previous value to avoid noise issues
                rmin(k) = sqrt(rmin(k)*r2)    
             end if

          end do
      end do
999   continue
      if (docache) then
         do l = 1, 3*N
            cache(l,1) = xx(l)
         end do
         do l = 1, 3*N
            cache(l,2) = vx(l)
         end do
         tcache = time
      end if
      RETURN
      END


c Implement effects of close encounter, returning whether we have merged or not
      LOGICAL FUNCTION implementFlyby(k,time,newreg,xx,vx,tmintau)
      INCLUDE "ARCCOM2e2.CH"
      include "rmins.h"
      include "orbel.h"
      LOGICAL newreg
      REAL*8 xx(*), vx(*), time
      INTEGER K

      real*8 xcom(3), vcom(3), dv1, dv2, dum,vkep,
     :       q1(3),q2(3),v1(3),v2(3),vsq,
     :       DELTAT,TMAX,stepr,soft,cmet(3),tolerance,mvec(3),
     :       m1,m2,m3,mt,mu,mu_out,mup,mmin
      integer idum,Nrun,IWR,Ixc,Nbh,i,j,l,i1,i2,i3,i1m,i2m,i3m,chkmul
      CHARACTER*22 outfile

c Import orbital elements etc. from fitting of r_min
      real*8 r,meanmot,rminout,Ean,truean,eccprime,
     :       ecc2prime,tmintau

      LOGICAL debug

      integer kpair(9)
      data kpair /1,2,3,1,3,2,2,3,1/

      debug = .false.
      implementFlyby = .false.

c Determine indices of colliding stars (i1,i2) and companion (i3)
      i1 = kpair(3*k-2)
      i2 = kpair(3*k-1)
      i3 = kpair(3*k)
      i1m = 3*(i1-1)
      i2m = 3*(i2-1)
      i3m = 3*(i3-1)
      open (27,FILE="debugFlybys.out",access="append")

      write (27,*) 'Implementing flyby', i1, i2, i3

c Neglect encounters that are not with the black hole
      if (i3.eq.3) then
         write (*,*) 'Unexpected item in bagging area!'
         return
      end if

c Obtain masses
      m1 = m(i1)
      m2 = m(i2)
      m3 = m(i3)
      mt = m1+m2
      mu = m1*m2/mt
      mu_out = mt*m3/(mt+m3)
      mvec(1) = m1
      mvec(2) = m2
      mvec(3) = m3

      mmin = min(m1,m2)

      if (debug) then
         write (21,*) 0,(xx(j),j=1,9),(vx(j),j=1,9)
      end if

      if (debug) then
         write (21,*) 1, (rvec(j),j=1,3),(vvec(j),j=1,3)
      end if

c Obtain velocity and position of COM (constant) and relative for the two stars pre-enc
      vsq = 0.d0
      do j = 1, 3
         xcom(j) = (m1*xx(i1m+j) + m2*xx(i2m+j))/mt
         q1(j) = cache(i1m+j,1) - xcom(j)
         q2(j) = cache(i2m+j,1) - xcom(j)
         vcom(j) = (m1*vx(i1m+j) + m2*vx(i2m+j))/mt
         v1(j) = cache(i1m+j,2) - vcom(j)
         v2(j) = cache(i2m+j,2) - vcom(j)
         vsq = vsq + (v1(j)-v2(j))**2
      end do

c XXX TODO perhaps re-think order of calculating these things
      h = sqrt(hvec(1)**2 + hvec(2)**2 + hvec(3)**2)

c Obtain pre-encounter orbital elements
c Now cached in orbel.h
c      call xvToOrbEl(rvec, vvec, hvec, mt,
c     :          semi, ecc, inc, longasc, longperi, tmintau)

      if (debug) then
         write (21,*) 2,semi,ecc,inc,longasc,longperi,tmintau,h
      end if

c Calculate r_min (q)
      q = semi*(1.d0-ecc)

c Obtain encounter parameters
      call defineEncounter(q, k, m1, m2, alphaJ, dm, dspe)
      write (27,*) 'Encounter:q m1 m2 aJ dm dspe',q,m1,m2,alphaJ,dm,dspe
      if (debug) then
         write (21,*) 3,q,m1,m2,alphaJ,dm,dspe
      end if

c Test for complete mass transfer (merger)
      if (dm.ge.mmin) then
         if (debug) write (21,*) 4,'Merger'
         implementFlyby = .true.
         return
      end if

c Identify mass losing star as lower mass of the two and calculate new masses
      if (m1.lt.m2) then
         m1p = m1 - dm
         m2p = m2 + dm
      else
         m1p = m1 + dm
         m2p = m2 - dm
      end if

c mu prime and alpha_mu
      mup = m1p*m2p/mt
      alphamu = mup/mu

c Obtain new orbital elements
      hprime = alphaJ/alphamu*h
      semiprime = semi/(1.d0 - 2.d0*dspe*semi/mt)
      ecc2prime = 1.d0-hprime**2/(mt*semiprime)
      if (ecc2prime.gt.0.d0) then
         eccprime = dsqrt(ecc2prime)
      else
         eccprime = ecc / 2.d0    ! Shrink but keep finite for calculational convenience
      end if

      if (debug) then
         write (21,*) 4,mup, alphamu,hprime,semiprime,eccprime
      end if

c Check for mass transfer in a substantially circular orbit (implies merger)
      if (ecc.lt.1.e-3.and.dm.gt.mmin/1.d1) then
         if (debug) write (21,*) 'Mass transfer in circular orbit'
         implementFlyby = .true.
         return
      end if


c Convert to XYZUVW
      call OrbElToxv(semiprime,eccprime,inc,longasc,longperi,tmintau,mt,
     :                    rvec,vvec)
      if (debug) then
         write (21,*) 5,(rvec(j),j=1,3),(vvec(j),j=1,3)
      end if


      close(27)

c Calculate the new positions and velocities
      do j = 1, 3
         xx(i1m+j) = xcom(j) + m2p/mt * rvec(j)
         xx(i2m+j) = xcom(j) - m1p/mt * rvec(j)
         vx(i1m+j) = vcom(j) + m2p/mt * vvec(j)
         vx(i2m+j) = vcom(j) - m1p/mt * vvec(j)
      end do
      m(i1) = m1p
      m(i2) = m2p

      if (debug) then
         write (21,*) 6,(xx(j),j=1,9),(vx(j),j=1,9)
      end if

      newreg = .true.

      RETURN
      END

c This routine implements all of the physics from the SPH calculations
      SUBROUTINE defineEncounter(q, kp, m1, m2, alphaJ, dm, dspe)
      IMPLICIT NONE
      include "ARCparams.h"
      include "rmins.h"

      integer kp, nd
      REAL*8 q, m1, m2, alphaJ, dm, dspe, q2, qm
      REAL*8 mmin, mmax, zz(4), dlq

      mmin = min(m1,m2)
      mmax = m1+m2-mmin

      q2 = q**2

c Decide how to define encounter
      if (q2.lt.rcoll(kp)) then
         write (*,*) 'Fatal error in defineEncounter: too close!', 
     :               q, kp, m1, m2, rcoll(kp), rtid(kp)
         stop
      end if

c Test for the case where we are off the beginning of the table
      if (q.lt.dexp(tidDat(1,1,kp))) then
         write (*,*) 'Warning: defineEncounter: off beginning of table',
     :      kp, q, dexp(tidDat(1,1,kp))
         dspe = -1.d0 * dexp(tidDat(1,2,kp))
         dm = dexp(tidDat(1,3,kp))
         alphaJ = 1.d0 - dexp(tidDat(1,4,kp))
         write (*,*) 'dspe dm aJ', dspe, dm, alphaJ
         return
      end if

c Test for the case where we are off the end of the table
      if (q2.gt.rtid(kp)) then
         write (*,*) 'Warning: defineEncounter: too far!', 
     :               q, kp, m1, m2, rcoll(kp), rtid(kp)
         dspe = 0.d0
         dm = 0.d0
         alphaJ = 1.d0
         return
      end if

c Extend the table with a -6 power law & exp cutoff (Fabian,Pringle&Rees)
      nd = ndx(kp)
c      write (*,*) kp, ndx(kp)
c      qmax = dexp(tidDat(nd,1,kp))
      qm = qmax(kp)
      if (q.gt.qm) then
         dlq = -6.d0*log(q/qm) - depl(kp)*(q**1.5-qm**1.5)
         dspe = -dexp(tidDat(nd,2,kp) + dlq)
         dm = 0.d0
         alphaJ = 1.d0 - dexp(tidDat(nd,4,kp) + dlq)
         write (27,*) 'Extension: kp nd qmax dlq q ', kp, nd, qm, dlq, q
         return
      end if

c Interpolate in tables
      call splint(tidDat(1,1,kp), tidSplin(1,1,kp), ndx(kp),
     :            4, log(q), zz)
      dspe = -1.d0 * dexp(zz(2))
      dm = dexp(zz(3))
      alphaJ = 1.d0 - dexp(zz(4))
      return
      

c 0.6 Mo results
      if (mmin.lt.0.601) then
         alphaJ = 1.d0
         dm = 0.d0
         dspe = -.02/(q/2.d0)**6
      else
c 0.8 Mo results
         alphaJ = 1.d0
         dm = 0.d0
         dspe = -.02/(q/2.d0)**6
      end if

      RETURN
      END


c Implement tidal capture and work out initial conditions to continue
      SUBROUTINE implementTidalCapt(k,time)
c      IMPLICIT NONE
      INCLUDE "ARCCOM2e2.CH"
      include "rmins.h"

      real*8 xcom(3), vcom(3), xx(9), vx(9), dv1, dv2, dum,vkep,
     :       DELTAT,TMAX,stepr,soft,cmet(3),tolerance,mvec(3),
     :       m1,m2,m3,mt,mu
      integer idum,Nrun,IWR,Ixc,Nbh,i,j,l,i1,i2,i3,i1m,i2m,i3m,chkmul
      CHARACTER*22 outfile
      integer kpair(9)
      data kpair /1,2,3,1,3,2,2,3,1/

c Determine indices of colliding stars (i1,i2) and companion (i3)
      i1 = kpair(3*k-2)
      i2 = kpair(3*k-1)
      i3 = kpair(3*k)
      i1m = 3*(i1-1)
      i2m = 3*(i2-1)
      i3m = 3*(i3-1)

c Obtain masses
      m1 = m(i1)
      m2 = m(i2)
      m3 = m(i3)
      mt = m1+m2
      mu = mt*m3/(mt+m3)
      mvec(1) = m1
      mvec(2) = m2
      mvec(3) = m3

c Obtain velocity and position of COM of post-capture binary
      do j = 1, 3
         xcom(j) = (m1*cache(i1m+j,1) + m2*cache(i2m+j,1))/mt
         vcom(j) = (m1*cache(i1m+j,2) + m2*cache(i2m+j,2))/mt
      end do

c Calculate r_min
      r = 0.d0
      do j = 1, 3
         r = r + (cache(i1m+j,1)-cache(i2m+j,1))**2
      end do
      r = sqrt(r)
      vkep = sqrt(0.5d0*mt/r)

c New positions
      do j = 1, 3
         xx(j)   = 2.d0*cache(i1m+j,1)-xcom(j)
         xx(3+j) = 2.d0*cache(i2m+j,1)-xcom(j)
         xx(6+j) = cache(i3m+j,1)
      end do

c Velocity differences from COM (this preserves the orbital plane)
      dv1 = 0.d0
      dv2 = 0.d0
      do j = 1, 3
         dv1 = dv1 + (cache(i1m+j,2)-vcom(j))**2
         dv2 = dv2 + (cache(i2m+j,2)-vcom(j))**2
      end do
      dv1 = sqrt(dv1)
      dv2 = sqrt(dv2)

c New velocities
      do j = 1, 3
         vx(j) = vcom(j) + (cache(i1m+j,2)-vcom(j))/dv1 * vkep * m2/mt
         vx(3+j) = vcom(j) + (cache(i2m+j,2)-vcom(j))/dv2 * vkep * m1/mt
         vx(6+j) = cache(i3m+j,2)
      end do

c Open cache of configuration options and read them back in
      open (81,file="crap.out")
      read (81,*) irun,IWR,idum,DELTAT,CHKMUL,TMAX,stepr,soft,cmet,
     &                dum,outfile,Ixc,Nbh ,dum,dum,dum,tolerance
      close(81)

      open (80, file="tidCap.out",access="append")
      write (80,*) irun, time, 2*r
      close(80)

c More detail in time after restart
      chkmul = chkmul / 10

c Write out configuration file to continue from
      write(outfile,'("cont.", i6.6, ".out")'), irun
      open (80,file="cnf.AR.cont",access="append")
      write (80,*) iwr, 3, deltat, chkmul, tmax,stepr,soft,cmet,
     &      Clight,outfile,Ixc,Nbh ,spin,tolerance
      do i = 1, 3
         L=3*(I-1)
         write(80,*)Mvec(I),(XX(L+K),K=1,3),(VX(L+K),K=1,3)
      end do
      close(80)
      return
      end
         

c 

c Implement collision and look at properties of resulting two-body system
c ictype=2 => not a collision but tidal capture (still, 2-body output may be useful)
      SUBROUTINE implementCollision(k,ictype)
      INCLUDE "ARCCOM2e2.CH"
      include "rmins.h"
      character*20 outfile
      real*8 xx(6),vx(6),dx(3),dv(3),m1,m2,m3,mt,mtt,mu,
     :       semi,ecc2,hvec(3),h,h2,cmet(3),pi,ecc
      integer kpair(9)
      data kpair /1,2,3,1,3,2,2,3,1/
      data pi /3.14159265358979324d0/

c Determine indices of colliding stars (i1,i2) and companion (i3)
      i1 = kpair(3*k-2)
      i2 = kpair(3*k-1)
      i3 = kpair(3*k)
      i1m = 3*(i1-1)
      i2m = 3*(i2-1)
      i3m = 3*(i3-1)

c Obtain masses
      m1 = m(i1)
      m2 = m(i2)
      m3 = m(i3)
      mt = m1+m2
      mtt = m3+m3
      mu = mt*m3/mtt

c Momentum and COM are conserved so we preserve the COM frame.

c Copy companion star
      do j = 1, 3
         xx(3+j) = cache(i3m+j,1)
         vx(3+j) = cache(i3m+j,2)
      end do

c  Create combined star in other place and obtain orbital elements
      v2 = 0.d0
      r2 = 0.d0
      do j = 1, 3
         xx(j) = (m1*cache(i1m+j,1) + m2*cache(i2m+j,1))/mt
         dx(j) = xx(j)-xx(j+3)
         x2 = x2 + dx(j)**2
         vx(j) = (m1*cache(i1m+j,2) + m2*cache(i2m+j,2))/mt
         dv(j) = vx(j)-vx(j+3)
         v2 = v2 + dv(j)**2
      end do

c Obtain orbital elements
      r = sqrt(x2)
      energy = .5d0*mu*v2 - m3*mt/r
      write (*,*) 'Energy=', energy, mu, v2, m3, mt, r
      semi = -(mt*m3)/(2.d0*energy)
      hvec(1) = dx(2)*dv(3)-dx(3)*dv(2)
      hvec(2) = dx(3)*dv(1)-dx(1)*dv(3)
      hvec(3) = dx(1)*dv(2)-dx(2)*dv(1)
      h2 = hvec(1)**2 + hvec(2)**2 + hvec(3)**2
      ecc2 = 1-h2/(mtt*semi)
      ecc = sqrt(ecc2)

c Write orbital elements out and return
      open (81,file="coll.out",access="append")
      write (81,99003) irun, ictype, i1, i2, i3, semi, ecc, 
     :                 semi*(1-ecc), r, (hvec(j),j=1,3), h2, mu
     :                  
      close(81)
99003 FORMAT(i6, 4(1x,i1), 1p, 9(1x,e15.3))

c XXX TODO REMOVE THIS lazy check of two-body orbit
      open (81,file="crap.out")
      read (81,*) irun,IWR,idum,DELTAT,CHKMUL,TMAX,stepr,soft,cmet,
     &                dum,outfile,Ixc,Nbh ,dum,dum,dum,tolerance
      close(81)

c Write out configuration file to continue from
      write(outfile,'("cont.", i6.6, ".out")'), irun
      open (80,file="cnf.AR.2body",access="append")
      tmax = 20.d0*pi*sqrt(semi**3/mtt)
      deltat = tmax/1.d3
      write (80,*) irun,iwr, 2, deltat, 10, tmax,stepr,soft,cmet,
     &      Clight,outfile,Ixc,Nbh ,spin,tolerance
      write(80,*)mt,(XX(j),j=1,3),(VX(j),j=1,3)
      write(80,*)m3,(XX(j),j=4,6),(VX(j),j=4,6)
      close(80)
      
      RETURN
      END
      


c Write out details of global closest encounters
      SUBROUTINE overallRmin
      INCLUDE "ARCCOM2e2.CH"
      include "rmins.h"
      CHARACTER*40 fn
      INTEGER i,j,k,l,kpair(9)

      data kpair /1,2,3,1,3,2,2,3,1/
      
      open (81, file="rmin.global.out", access="append")
      write (81,99002) (sqrt(rminglobal(k)),k = 1, N*(N-1)/2)
      close(81)
99002 FORMAT(1p,40(e10.3,1x))
      do k = 1, N*(N-1)/2
         i = kpair(3*k-2)
         j = kpair(3*k-1)
         WRITE(fn,'("rmin.global.", i1, ".out")')  k
         OPEN (81,file=fn,access="append")
         write(81,99000) irun,tminglobal(k),i,j,sqrt(rminglobal(k)),
     :                (cacheglobal(l,1,k),l=1,3*N),
     :                (cacheglobal(l,2,k),l=1,3*N),
     :                (mcache(l),l=1,N)
99000 FORMAT(i6,1x,1p,e13.5,0p,2(1x,i2),1x,1p,75(e13.5,1x))
         CLOSE(81)
      end do

      RETURN
      END

c Initialise data for dealing with stellar collisions
      SUBROUTINE initRmin
      INCLUDE "ARCCOM2e2.CH"
      INTEGER N
      character*40 fn
c
      INTEGER i,j,k,np,nskip
      real*8 RM1, RC1
c
      include "rmins.h"
c
c Set up large initial values for rmin
      k = 0
      do i = 1, N-1
         do j = i+1, N
            k = k + 1
            rmin(k) = 1.d20
            rminglobal(k) = 1.d20
            tminglobal(k) = -1.d0
            printed(k) = .FALSE.
         end do
      end do
      np = k
c
c Read in data to deal with interactions between pairs of bodies
      k = 0
      do i = 1, N-1
         do j = i+1, N
            k = k + 1
            write(fn,'("tidDat.", i1, ".", i1, ".dat")') i, j
            open(81,err=10,file=fn,status="old")
            read (81,*) ! Header line
            read (81,*) rtid(k), rcoll(k), ndx(k), depl(k)
            if (rtid(k).gt.rcoll(k).and.ndx(k).lt.1) then
               write (*,*) 'Fatal error reading tidal parameters'
               write (*,*) k, rtid(k), rcoll(k), depl(k)
               STOP
            end if
            rtid(k) = rtid(k)**2
            rcoll(k) = rcoll(k)**2
            nskip = 0
            if (ndx(k).eq.0) goto 20
            do kk = 1, ndx(k)
               read (81,*) q, de, dm, alphaJ
               if (q.gt.0.d0.and.de.lt.0.d0.and.dm.gt.0.d0.and.alphaJ.lt.1.d0) then
                  tidDat(kk,1,k) = log(q)
                  tidDat(kk,2,k) = log(-de)
                  tidDat(kk,3,k) = log(dm)
                  tidDat(kk,4,k) = log(1-alphaJ)
               else
                  nskip = nskip + 1
                  write (*,*) 'Stars ', i, j, ' pair ', k,
     :               ' skipped tidal datum', q, de, dm, alphaJ
               end if
            end do
            ndx(k) = ndx(k) - nskip
            close(81)

c Store furthest out point in the table
            qmax(k) = dexp(tidDat(ndx(k),1,k))

c Use the data that we have loaded to build splines
            call spline(tidDat(1,1,k), ndx(k), 4, tidSplin(1,1,k))
            goto 20

c  No file found -- non-interacting particles
10          ndx(k) = 0
            rcoll(k) = 0.d0
            rtid(k) = 0.d0
            qmax(k) = 0.d0
20          continue
         end do
      end do

      
c
      RETURN
      END

c Work out the correct values to use for r_min based on the orbital properties
      SUBROUTINE fitrmin(i,j,k,XX,VX,rminout,eccout, tmintau)
      INCLUDE "ARCCOM2e2.CH"
      INCLUDE "rmins.h"
      INCLUDE "orbel.h"
      INTEGER i,j,k
      REAL*8 XX(*),VX(*)

      REAL*8 R2,V2,h2,xk,xj,vk,vj,mt,r,
     :       meanmot,Ean,rminout,eccout,
     :          h
      INTEGER l
      

c Calculate energy and angular momentum
      R2 = 0.d0
      V2 = 0.d0
      mt = m(i) + m(j)
      do l = 1, 3
         xk = XX(3*i-3+l)
         xj = XX(3*j-3+l)
         vk = VX(3*i-3+l)
         vj = VX(3*j-3+l)
         rvec(l) = xk-xj
         vvec(l) = vk-vj
         R2 = R2 + (xk-xj)**2
         V2 = V2 + (vk-vj)**2
      end do
      hvec(1) = rvec(2)*vvec(3)-rvec(3)*vvec(2)
      hvec(2) = rvec(3)*vvec(1)-rvec(1)*vvec(3)
      hvec(3) = rvec(1)*vvec(2)-rvec(2)*vvec(1)
      h2 = hvec(1)**2 + hvec(2)**2 + hvec(3)**2
      h = dsqrt(h2)

      call xvToOrbEl(rvec, vvec, hvec, mt,
     :          semi, ecc, inc, longasc, longperi, tmintau)

      rminout = semi*(1.d0-ecc)
      eccout = ecc

c      write (22,*) 'fitrmin', irun, i, j, k, 'r', (rvec(l),l=1,3),
c     :             'v', (vvec(l),l=1,3), 'h', (hvec(l),l=1,3),
c     :             r2,v2,mt, 'a e q', semi, ecc, rminout,
c     :             'dt', tmintau
      return
      end


c Convert orbital elements to Cartesians using John Chambers' routine
      SUBROUTINE OrbElToxv(semi,ecc,inc,longasc,longperi,
     :                    tmintau,mt,rvec,vvec)
      IMPLICIT NONE
      REAL*8 semi,ecc,inc,longasc,longperi,tmintau,mt,rvec(3),vvec(3),
     :       x,y,z,u,v,w,q,meanan,meanmot
      
      q = semi*(1.d0-ecc)
      meanmot = sqrt(abs(mt/semi**3))
      meanan = tmintau * meanmot

      call mco_el2x(mt,q,ecc,inc,longperi,longasc,meanan,
     :              rvec(1),rvec(2),rvec(3),vvec(1),vvec(2),vvec(3))

      return
      end

c Convert cartesians to orbital elements using John Chambers' routine
      SUBROUTINE xvtoorbel(rvec, vvec, hvec, mt,
     :          semi, ecc, inc, longasc, longperi, tmintau)

      IMPLICIT NONE
      REAL*8 rvec(3),vvec(3),hvec(3),mt,semi,ecc,inc,longasc,longperi,
     :       truean,tmintau

      REAL*8 q,meanmot,meanan

      CALL mco_x2el(mt,rvec(1),rvec(2),rvec(3),vvec(1),vvec(2),vvec(3),
     :            q,ecc,inc,longperi,longasc,meanan)

      semi = q/(1.d0-ecc)
      meanmot = sqrt(abs(mt/semi**3))
      tmintau = meanan/meanmot
      
      return
      end

