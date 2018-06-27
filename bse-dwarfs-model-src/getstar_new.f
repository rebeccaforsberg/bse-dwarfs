c Obtain a mass from a number of imfs
      real*8 function getmass(imf, mmin, mmax)
      implicit none
      integer imf

      include "rand.h"

      real*8 X,imf1,imf2,imf3, mmin, mmax
      real*8 ran3

      getmass = 0.d0
      do while (getmass.lt.mmin .or. getmass.gt.mmax) 
         X = ran3(seed)
         if (imf.eq.0) then
            getmass = mmin
         else if (imf.eq.1) then
            getmass = imf1(X)
         else if (imf.eq.2) then
            getmass = imf2(X, mmin, 1.35d0)
         else
            getmass = imf3(X)
         end if
      end do

      return
      end

      real* 8 function dm91(m1)
      implicit none

      include "rand.h"

      real*8 gasdev, logp

      if (m1.gt.0.1 .and. m1.le.0.8)
         LOGP = 
      else if (m1.gt. 0.8 .and m1.le.1.6) then
*       Utilise Duquennoy & Mayor (1991) distribution of periods
         LOGP = 2.3D0*GASDEV() + 5D0
c Period in days

      else if (m1.gt.1.6 .and. m1.le.5) then

      end if   
      dm91 = 10.D0**LOGP
      return
      end

      function r10()
      implicit none
      real*8 r10

      include "rand.h"

      real*8 gasdev, logp

*       Utilise Raghavan (2010) distribution of periods
      LOGP = 2.28D0*GASDEV() + 5.03D0
c Period in days
      r10 = 10.D0**LOGP

      return
      end

      function opik(amin,amax)
      implicit none
      real*8 opik,amin,amax,ran3

      include "rand.h"

      opik = amin * dexp(ran3(seed) * log(amax/amin))
      return
      end




