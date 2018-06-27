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

      function dm91()
      implicit none
      real*8 dm91

      include "rand.h"

      real*8 gasdev, logp

*       Utilise Duquennoy & Mayor (1991) distribution of periods
      LOGP = 2.3D0*GASDEV() + 4.8D0
c Period in days
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




