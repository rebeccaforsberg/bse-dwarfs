c Power law mass ratio
 
c To do: Add periods
      SUBROUTINE masspowerlaw(m1,m2,tb)
      IMPLICIT NONE
      REAL*8 m1,m2,gammaq1,gammaq2,qb,X

      qb = 0.3d0

      include "rand.h"
      REAL*8 ran3

      if (m1.ge.0.8d0 .and. m1.lt.1.6d0) then
            gammaq2 = 0.3d0
            if (log10(tb).gt.0.d0 .and. log10(tb).le.6.d0) then
                  gammaqh = -0.5d0
            else if (log10(tb).gt.6.d0) then
                  gammaql = -1.1d0
            end if
c What to do in gap???
      else if (m1.ge.6.d0 .and. m1.lt.5.d0) then
            if (log10(tb).gt.0.d0 .and. log10(tb).le.2.d0) then
                  gammaql = -0.5d0
                  gammaqs =  0.2d0
            else if (log10(tb).gt.2.d0 .and. log10(tb).le.4.d0) then
                  gammaql = -0.9d0
                  gammaqs =  0.1d0
            else if (log10(tb).gt.4.d0 .and. log10(tb).le.6.d0) then
                  gammaql = -1.4d0
                  gammaqs = -0.5d0
            else if (log10(tb).gt.6.d0) then
                  gammaql = -2.d0
                  gammaqs = -1.d0
            end if

      else if (m1.ge.5.d0) then
            if (log10(tb).gt.0.d0 .and. log10(tb).le.2.d0) then
                  gammaql = -0.5d0
                  gammaqs =  0.1d0
            else if (log10(tb).gt.2.d0 .and. log10(tb).le.4.d0) then
                  gammaql = -1.7d0
                  gammaqs = -0.2d0
            else if (log10(tb).gt.4.d0 .and. log10(tb).le.6.d0) then
                  gammaql = -2.0d0
                  gammaqs = -1.2d0
            else if (log10(tb).gt.6.d0) then
                  gammaql = -2.0d0
                  gammaqs = -1.5d0
            end if
      end if

      X = ran3
      
      m2 = m1*

      RETURN
      END
