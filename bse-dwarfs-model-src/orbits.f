      PROGRAM testOrbEl
      IMPLICIT NONE
      INTEGER j
      REAL*8 a, m1, m2, dx, dv

      INCLUDE "rand.h"

      seed = -12345

      a = 1.d0
      m1 = .5d0
      m2 = .5d0


      open (10,file="out.circ")
      open (11,file="out.e.5")
      open (12,file="out.e.9")
      do j = 1, 1000
         call findSepRvel(a, 0.d0,  m1, m2, dx, dv)
         write (10,*) j, dx, abs(dv)
         call findSepRvel(a, 5.d-1, m1, m2, dx, dv)
         write (11,*) j, dx, abs(dv)
         call findSepRvel(a, 9.d-1, m1, m2, dx, dv)
         write (12,*) j, dx, abs(dv)
      end do
      close(10)
      close(11)
      close(12)

      END



