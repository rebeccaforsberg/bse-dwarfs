      PROGRAM sampleGalaxy
      IMPLICIT NONE
      REAL*8 Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal
      integer i,Nsamp
      logical galReject

      parameter(Nsamp=400000)

      include "rand.h"

      seed = -54321

      open (20,file="galaxy.out",status="new")
      write (20,99001)

c Draw 1e4 samples from the whole Galaxy
      do i = 1, Nsamp
         call galacticModel(Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal)
         write (20,99002) i, Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal
      end do

      close(20)

c Draw 1e4 samples from the segment of the Galaxy surveyed by APOGEE
      open (20,file="apogeeGalaxy.out",status="new")
      write (20,99001)
      do i = 1, Nsamp
         galReject = .true.
         do while (galReject)
            call galacticModel(Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal)
c Selection criteria to mimic APOGEE
            if (lgal.lt.24.d0 .or. lgal.gt.24.d1 .or. abs(bgal).gt.16.d0) then
               galReject = .true.
            else
               galReject = .false.
            end if
         end do
         write (20,99002) i, Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal
      end do
      close(20)
         

99001 format("#i        R        q        z        D        l        b        x        y")
99002 format(i6,8(1x,f8.4))

      end
