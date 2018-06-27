c Implement a simple Galactic model
      SUBROUTINE galacticModel(R,q,z,D,l,b,x,y)
      IMPLICIT NONE
      REAL*8 R,q,z,D,l,b,x,y

      include "rand.h"
      REAL*8 ran3

      REAL*8 xSun
      REAL*8 Rdisc
      REAL*8 zdisc

      REAL*8 pi,twopi, dtor
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
      dtor = twopi / 360.d0

      xSun = -8.d0
      Rdisc = 3.d0
      zdisc = 0.3d0


c Positions in the disc (Galactocentric)
      R = -Rdisc * log(ran3(seed))
      z = zdisc * log(ran3(seed))
      if (ran3(seed).gt.0.5d0) z =  -1.d0 * z
      q = twopi * ran3(seed)
      x = R * sin(q)
      y = R*cos(q)
      q = q/dtor

c Solar-centric co-ordinates (D,l,b)
      D = sqrt((x-xsun)**2 + y**2 + z**2)
      b = asin(z/D) / dtor
      l = atan2(y,x-xsun) / dtor
      if (l.lt.0.d0) l = l + 360.d0

      RETURN
      END

