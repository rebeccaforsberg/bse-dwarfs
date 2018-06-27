c Find separation on the sky and relative velocity along the line of sight for a binary
c semi should be in rsun and m1, m2 in solar masses
c Then xdiff will be in rsun and vdiff in km/s
      SUBROUTINE findSepRvels(semi, ecc, m1, m2, times, ntimes,
     :                        xdiffs, vdiffs)
      IMPLICIT NONE
      REAL*8 semi,ecc,m1,m2,times(*),xdiffs(*),vdiffs(*)
      INTEGER ntimes
      REAL*8 q,meanmot,meanan,inc,longasc,longperi,mt,rvec(3),vvec(3),
     :       vunit,pi,ran3,meanan0,tunits,tunity, tunitd
      INTEGER j

      PARAMETER (vunit=436.82476)
      PARAMETER (tunits=1592.1716)            ! Derived time unit in s (== dynamical timescale for Sun)
      PARAMETER (tunity=5.04539537e-5)        ! Derived time unit in yr
      PARAMETER (tunitd=0.018427344)          ! Derived time unit in days
      PARAMETER (pi=3.1415926536)

      INCLUDE "rand.h"

c Derived quantities
      mt = m1 + m2
      q = semi*(1.d0-ecc)

      meanmot = sqrt(abs(mt/semi**3))

c Randomly generated quantities: mean anomoly uniform between 0 and 2*pi
      meanan = ran3(seed) * 2.d0 * pi
      meanan0 = meanan
c Longitude of ascending node uniform in [0:2pi]
      longasc = ran3(seed) * 2.d0 * pi
c Longitude of pericentre uniform in [0:2pi]
      longperi = ran3(seed) * 2.d0 * pi
c Inclination between 0 and pi with p(i)=.5 sin(i)
      inc = acos(1.d0-2.d0*ran3(seed))

c Loop over observation times
      do j = 1, ntimes
c         meanan = meanan0 + meanmot * (times(j)-times(1)) / tunity
         meanan = meanan0 + meanmot * (times(j)-times(1)) / tunitd
c Convert orbital elements to Cartesians using John Chambers' routine
         call mco_el2x(mt,q,ecc,inc,longperi,longasc,meanan,
     :              rvec(1),rvec(2),rvec(3),vvec(1),vvec(2),vvec(3))

c Sky is the xy plane so take separation in that
         xdiffs(j) = sqrt(rvec(1)**2 + rvec(2)**2)

c Radial velocity is along the z direction
         vdiffs(j) = vvec(3) * vunit
      end do


c      write (*,*) meanan, longperi, longasc, inc, mt
c      write (*,*) (rvec(j),j=1,3),(vvec(j),j=1,3), vunit
c      stop

      return
      end

c Find separation on the sky and relative velocity along the line of sight for a binary
c semi should be in rsun and m1, m2 in solar masses
c Then xdiff will be in rsun and vdiff in km/s
      SUBROUTINE findSepRvel(semi, ecc, m1, m2,
     :                        xdiff, vdiff)
      IMPLICIT NONE
      REAL*8 semi,ecc,m1,m2,xdiff,vdiff
      REAL*8 q,meanmot,meanan,inc,longasc,longperi,mt,rvec(3),vvec(3),
     :       vunit,pi,ran3
      INTEGER j

      PARAMETER (vunit=436.82476)
      PARAMETER (pi=3.1415926536)

      INCLUDE "rand.h"
c     write (*,*)'Here we are'
c Derived quantities
      mt = m1 + m2
      q = semi*(1.d0-ecc)
c      write (*,*)'Here we are now: mt, semi:', mt, semi
      meanmot = sqrt(abs(mt/semi**3))
c      write (*,*)'Here we are Last'
c Randomly generated quantities: mean anomoly uniform between 0 and 2*pi
      meanan = ran3(seed) * 2.d0 * pi
c Longitude of ascending node uniform in [0:2pi]
      longasc = ran3(seed) * 2.d0 * pi
c Longitude of pericentre uniform in [0:2pi]
      longperi = ran3(seed) * 2.d0 * pi
c Inclination between 0 and pi with p(i)=.5 sin(i)
      inc = acos(1.d0-2.d0*ran3(seed))
      
c Convert orbital elements to Cartesians using John Chambers' routine
      call mco_el2x(mt,q,ecc,inc,longperi,longasc,meanan,
     :              rvec(1),rvec(2),rvec(3),vvec(1),vvec(2),vvec(3))

c Sky is the xy plane so take separation in that
      xdiff = sqrt(rvec(1)**2 + rvec(2)**2)

c Radial velocity is along the z direction
      vdiff = vvec(3) * vunit

c      write (*,*) meanan, longperi, longasc, inc, mt
c      write (*,*) (rvec(j),j=1,3),(vvec(j),j=1,3), vunit
c      stop

      return
      end

c-------------------------------------------------------------------
c
***********************************************************************
c                    ORBEL_FGET.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_fget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
*           Cel. Mech. ".  Quartic convergence from Danby's book.
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: 2/26/93 hfl
*     Modified by JEC
***********************************************************************

      real*8 function orbel_fget(e,capn)

      include 'swift.inc'

c...  Inputs Only: 
      real*8 e,capn

c...  Internals:
      integer i,IMAX
      real*8 tmp,x,shx,chx
      real*8 esh,ech,f,fp,fpp,fppp,dx
      PARAMETER (IMAX = 10)

c----
c...  Executable code 

c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. 

c  begin with a guess proposed by Danby
      if( capn .lt. 0.d0) then
         tmp = -2.d0*capn/e + 1.8d0
         x = -log(tmp)
      else
         tmp = +2.d0*capn/e + 1.8d0
         x = log( tmp)
      endif

      orbel_fget = x

      do i = 1,IMAX
        shx = sinh(x)
        chx = cosh(x)
        esh = e*shx
        ech = e*chx
        f = esh - x - capn
c        write(6,*) 'i,x,f : ',i,x,f
        fp = ech - 1.d0  
        fpp = esh 
        fppp = ech 
        dx = -f/fp
        dx = -f/(fp + dx*fpp/2.d0)
        dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
        orbel_fget = x + dx
c   If we have converged here there's no point in going on
        if(abs(dx) .le. TINY) RETURN
        x = orbel_fget
      enddo

      write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
      return
      end   ! orbel_fget
c------------------------------------------------------------------
c
***********************************************************************
c                    ORBEL_FHYBRID.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                           n ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
*               For larger N, uses FGET
*     REMARKS: 
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26,1992.
*     REVISIONS: 
*     REVISIONS: 2/26/93 hfl
***********************************************************************

      real*8 function orbel_fhybrid(e,n)

      include 'swift.inc'

c...  Inputs Only: 
      real*8 e,n

c...  Internals:
      real*8 abn
        real*8 orbel_flon,orbel_fget

c----
c...  Executable code 

      abn = n
      if(n.lt.0.d0) abn = -abn

      if(abn .lt. 0.636d0*e -0.6d0) then
        orbel_fhybrid = orbel_flon(e,n)
      else 
        orbel_fhybrid = orbel_fget(e,n)
      endif   

      return
      end  ! orbel_fhybrid
c-------------------------------------------------------------------
c
***********************************************************************
c                    ORBEL_FLON.F
***********************************************************************
*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
*
*             Input:
*                           e ==> eccentricity anomaly. (real scalar)
*                        capn ==> hyperbola mean anomaly. (real scalar)
*             Returns:
*                  orbel_flon ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: Uses power series for N in terms of F and Newton,s method
*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 26, 1992.
*     REVISIONS: 
***********************************************************************

      real*8 function orbel_flon(e,capn)

      include 'swift.inc'

c...  Inputs Only: 
      real*8 e,capn

c...  Internals:
      integer iflag,i,IMAX
      real*8 a,b,sq,biga,bigb
      real*8 x,x2
      real*8 f,fp,dx
      real*8 diff
      real*8 a0,a1,a3,a5,a7,a9,a11
      real*8 b1,b3,b5,b7,b9,b11
      PARAMETER (IMAX = 10)
      PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
      PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
      PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
      PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

c----
c...  Executable code 


c Function to solve "Kepler's eqn" for F (here called
c x) for given e and CAPN. Only good for smallish CAPN 

      iflag = 0
      if( capn .lt. 0.d0) then
         iflag = 1
         capn = -capn
      endif

      a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
      a0 = -6227020800.d0*capn/e
      b1 = a1

c  Set iflag nonzero if capn < 0., in which case solve for -capn
c  and change the sign of the final answer for F.
c  Begin with a reasonable guess based on solving the cubic for small F


      a = 6.d0*(e-1.d0)/e
      b = -6.d0*capn/e
      sq = sqrt(0.25*b*b +a*a*a/27.d0)
      biga = (-0.5*b + sq)**0.3333333333333333d0
      bigb = -(+0.5*b + sq)**0.3333333333333333d0
      x = biga + bigb
c      write(6,*) 'cubic = ',x**3 +a*x +b
      orbel_flon = x
c If capn is tiny (or zero) no need to go further than cubic even for
c e =1.
      if( capn .lt. TINY) go to 100

      do i = 1,IMAX
        x2 = x*x
        f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
        fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
        dx = -f/fp
c        write(6,*) 'i,dx,x,f : '
c        write(6,432) i,dx,x,f
432     format(1x,i3,3(2x,1p1e22.15))
        orbel_flon = x + dx
c   If we have converged here there's no point in going on
        if(abs(dx) .le. TINY) go to 100
        x = orbel_flon
      enddo

c Abnormal return here - we've gone thru the loop 
c IMAX times without convergence
      if(iflag .eq. 1) then
         orbel_flon = -orbel_flon
         capn = -capn
      endif
      write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
        diff = e*sinh(orbel_flon) - orbel_flon - capn
        write(6,*) 'N, F, ecc*sinh(F) - F - N : '
        write(6,*) capn,orbel_flon,diff
      return

c  Normal return here, but check if capn was originally negative
100   if(iflag .eq. 1) then
         orbel_flon = -orbel_flon
         capn = -capn
      endif

      return
      end     ! orbel_flon
c
***********************************************************************
c                    ORBEL_ZGET.F
***********************************************************************
*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
*          given Q (Fitz. notation.)
*
*             Input:
*                           q ==>  parabola mean anomaly. (real scalar)
*             Returns:
*                  orbel_zget ==>  eccentric anomaly. (real scalar)
*
*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
*     REMARKS: For a parabola we can solve analytically.
*     AUTHOR: M. Duncan 
*     DATE WRITTEN: May 11, 1992.
*     REVISIONS: May 27 - corrected it for negative Q and use power
*            series for small Q.
***********************************************************************

      real*8 function orbel_zget(q)

      include 'swift.inc'

c...  Inputs Only: 
      real*8 q

c...  Internals:
      integer iflag
      real*8 x,tmp

c----
c...  Executable code 

      iflag = 0
      if(q.lt.0.d0) then
        iflag = 1
        q = -q
      endif

      if (q.lt.1.d-3) then
         orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
      else
         x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
         tmp = x**(1.d0/3.d0)
         orbel_zget = tmp - 1.d0/tmp
      endif

      if(iflag .eq.1) then
           orbel_zget = -orbel_zget
         q = -q
      endif
      
      return
      end    ! orbel_zget
c----------------------------------------------------------------------


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_EL2X.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates Cartesian coordinates and velocities given Keplerian orbital
c elements (for elliptical, parabolic or hyperbolic orbits).
c
c Based on a routine from Levison and Duncan's SWIFT integrator.
c
c  gm = grav const * (central + secondary mass)
c  q = perihelion distance
c  e = eccentricity
c  i = inclination                 )
c  p = longitude of perihelion !!! )   in
c  n = longitude of ascending node ) radians
c  l = mean anomaly                )
c
c  x,y,z = Cartesian positions  ( units the same as a )
c  u,v,w =     "     velocities ( units the same as sqrt(gm/a) )
c
c------------------------------------------------------------------------------
c
      subroutine mco_el2x (gm,q,e,i,p,n,l,x,y,z,u,v,w)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      real*8 gm,q,e,i,p,n,l,x,y,z,u,v,w
c
c Local
      real*8 g,a,ci,si,cn,sn,cg,sg,ce,se,romes,temp
      real*8 z1,z2,z3,z4,d11,d12,d13,d21,d22,d23
      real*8 mco_kep, orbel_fhybrid, orbel_zget
c
c------------------------------------------------------------------------------
c
c Change from longitude of perihelion to argument of perihelion
      g = p - n
c
c Rotation factors
      si = sin(i)
      ci = cos(i)
      sg = sin(g)
      cg = cos(g)
      sn = sin(n)
      cn = cos(n)
      z1 = cg * cn
      z2 = cg * sn
      z3 = sg * cn
      z4 = sg * sn
      d11 =  z1 - z4*ci
      d12 =  z2 + z3*ci
      d13 = sg * si
      d21 = -z3 - z2*ci
      d22 = -z4 + z1*ci
      d23 = cg * si
c
c Semi-major axis
      a = q / (1.d0 - e)
c
c Ellipse
      if (e.lt.1.d0) then
        romes = sqrt(1.d0 - e*e)
        temp = mco_kep (e,l)
        se = sin(temp)
        ce = cos(temp)
        z1 = a * (ce - e)
        z2 = a * romes * se
        temp = sqrt(gm/a) / (1.d0 - e*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
      else
c Parabola
        if (e.eq.1.d0) then
          ce = orbel_zget(l)
          z1 = q * (1.d0 - ce*ce)
          z2 = 2.d0 * q * ce
          z4 = sqrt(2.d0*gm/q) / (1.d0 + ce*ce)
          z3 = -ce * z4
        else
c Hyperbola
          romes = sqrt(e*e - 1.d0)
          temp = orbel_fhybrid(e,l)
          se = sinh(temp)
          ce = cosh(temp)
          z1 = a * (ce - e)
          z2 = -a * romes * se
          temp = sqrt(gm/abs(a)) / (e*ce - 1.d0)
          z3 = -se * temp
          z4 = romes * ce * temp
        end if
      endif
c
      x = d11 * z1  +  d21 * z2
      y = d12 * z1  +  d22 * z2
      z = d13 * z1  +  d23 * z2
      u = d11 * z3  +  d21 * z4
      v = d12 * z3  +  d22 * z4
      w = d13 * z3  +  d23 * z4
c
c------------------------------------------------------------------------------
c
      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_X2EL.FOR    (ErikSoft  23 January 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Calculates Keplerian orbital elements given relative coordinates and
c velocities, and GM = G times the sum of the masses.
c
c The elements are: q = perihelion distance
c                   e = eccentricity
c                   i = inclination
c                   p = longitude of perihelion (NOT argument of perihelion!!)
c                   n = longitude of ascending node
c                   l = mean anomaly (or mean longitude if e < 1.e-8)
c
c------------------------------------------------------------------------------
c
      subroutine mco_x2el (gm,x,y,z,u,v,w,q,e,i,p,n,l)
c
      implicit none
      include 'mercury.inc'
c
c Input/Output
      real*8 gm,q,e,i,p,n,l,x,y,z,u,v,w
c
c Local
      real*8 hx,hy,hz,h2,h,v2,r,rv,s,true
      real*8 ci,to,temp,tmp2,bige,f,cf,ce
c
c------------------------------------------------------------------------------
c
      hx = y * w  -  z * v
      hy = z * u  -  x * w
      hz = x * v  -  y * u
      h2 = hx*hx + hy*hy + hz*hz
      v2 = u * u  +  v * v  +  w * w
      rv = x * u  +  y * v  +  z * w
      r = sqrt(x*x + y*y + z*z)
      h = sqrt(h2)
      s = h2 / gm
c
c Inclination and node
      ci = hz / h
      if (abs(ci).lt.1) then
        i = acos (ci)
        n = atan2 (hx,-hy)
        if (n.lt.0) n = n + TWOPI
      else
        if (ci.gt.0) i = 0.d0
        if (ci.lt.0) i = PI
        n = 0.d0
      end if
c
c Eccentricity and perihelion distance
      temp = 1.d0  +  s * (v2 / gm  -  2.d0 / r)
      if (temp.le.0) then
        e = 0.d0
      else
        e = sqrt (temp)
      end if
      q = s / (1.d0 + e)
c
c True longitude
      if (hy.ne.0) then
        to = -hx/hy
        temp = (1.d0 - ci) * to
        tmp2 = to * to
        true = atan2((y*(1.d0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
      else
        true = atan2(y * ci, x)
      end if
      if (ci.lt.0) true = true + PI
c
      if (e.lt.3.d-8) then
        p = 0.d0
        l = true
      else
        ce = (v2*r - gm) / (e*gm)
c
c Mean anomaly for ellipse
        if (e.lt.1) then
          if (abs(ce).gt.1) ce = sign(1.d0,ce)
          bige = acos(ce)
          if (rv.lt.0) bige = TWOPI - bige
          l = bige - e*sin(bige)
        else
c
c Mean anomaly for hyperbola
          if (ce.lt.1) ce = 1.d0
          bige = log( ce + sqrt(ce*ce-1.d0) )
          if (rv.lt.0) bige = - bige
          l = e*sinh(bige) - bige
        end if
c
c Longitude of perihelion
        cf = (s - r) / (e*r)
        if (abs(cf).gt.1) cf = sign(1.d0,cf)
        f = acos(cf)
        if (rv.lt.0) f = TWOPI - f
        p = true - f
        p = mod (p + TWOPI + TWOPI, TWOPI)
      end if
c
      if (l.lt.0) l = l + TWOPI
      if (l.gt.TWOPI) l = mod (l, TWOPI)
c
c------------------------------------------------------------------------------
c
      return
      end


c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MCO_KEP.FOR    (ErikSoft  7 July 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Solves Kepler's equation for eccentricities less than one.
c Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
c
c  e = eccentricity
c  l = mean anomaly      (radians)
c  u = eccentric anomaly (   "   )
c
c------------------------------------------------------------------------------
c
      function mco_kep (e,oldl)
      implicit none
c
c Input/Outout
      real*8 oldl,e,mco_kep
c
c Local
      real*8 l,pi,twopi,piby2,u1,u2,ome,sign
      real*8 x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
      real*8 p,q,p2,ss,cc
      logical flag,big,bigg
c
c------------------------------------------------------------------------------
c
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
      piby2 = .5d0 * pi
c
c Reduce mean anomaly to lie in the range 0 < l < pi
      if (oldl.ge.0) then
        l = mod(oldl, twopi)
      else
        l = mod(oldl, twopi) + twopi
      end if
      sign = 1.d0
      if (l.gt.pi) then
        l = twopi - l
        sign = -1.d0
      end if
c
      ome = 1.d0 - e
c
      if (l.ge..45d0.or.e.lt..55d0) then
c
c Regions A,B or C in Nijenhuis
c -----------------------------
c
c Rough starting value for eccentric anomaly
        if (l.lt.ome) then
          u1 = ome
        else
          if (l.gt.(pi-1.d0-e)) then
            u1 = (l+e*pi)/(1.d0+e)
          else
            u1 = l + e
          end if
        end if
c
c Improved value using Halley's method
        flag = u1.gt.piby2
        if (flag) then
          x = pi - u1
        else
          x = u1
        end if
        x2 = x*x
        sn = x*(1.d0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.d0 + x2*(-.49815 + x2*.03805)
        if (flag) dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.d0 - e*dsn
        u2 = u1 - f0/(f1 - .5d0*f0*f2/f1)
      else
c
c Region D in Nijenhuis
c ---------------------
c
c Rough starting value for eccentric anomaly
        z1 = 4.d0*e + .5d0
        p = ome / z1
        q = .5d0 * l / z1
        p2 = p*p
        z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.d0*q / ( z2 + p + p2/z2 )
c
c Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3)
        u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 )
      end if
c
c Accurate value using 3rd-order version of Newton's method
c N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
c
c First get accurate values for u2 - sin(u2) and 1 - cos(u2)
      bigg = (u2.gt.piby2)
      if (bigg) then
        z3 = pi - u2
      else
        z3 = u2
      end if
c
      big = (z3.gt.(.5d0*piby2))
      if (big) then
        x = piby2 - z3
      else
        x = z3
      end if
c
      x2 = x*x
      ss = 1.d0
      cc = 1.d0
c
      ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. -
     %   x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
      cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. -
     %   x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. -
     %   x2/306.))))))))
c
      if (big) then
        z1 = cc + z3 - 1.d0
        z2 = ss + z3 + 1.d0 - piby2
      else
        z1 = ss
        z2 = cc
      end if
c
      if (bigg) then
        z1 = 2.d0*u2 + z1 - pi
        z2 = 2.d0 - z2
      end if
c
      f0 = l - u2*ome - e*z1
      f1 = ome + e*z2
      f2 = .5d0*e*(u2-z1)
      f3 = e/6.d0*(1.d0-z2)
      z1 = f0/f1
      z2 = f0/(f2*z1+f1)
      mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
c
c------------------------------------------------------------------------------
c
      return
      end
