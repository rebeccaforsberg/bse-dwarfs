c Load a metallicity distribution from disk
      subroutine loadzdist
      implicit none
      include "zdist.h"

      integer i

      open(35,file="zdist.dat",status="old")
      i = 0
      do
         i = i + 1
         read(35,*,end=100) zdist(i)
      end do
100   nzdist = i
      close(35)
      return
      end

c Obtain a metallicity from a metallicity distribution
      real*8 function getzcum(zmin, zmax, zfuzz)
      implicit none
      real*8 zmin,zmax,zfuzz
      include "zdist.h"
      include "rand.h"

      real*8 feh,X,ran3
      integer i

      feh = zmin - 1.d0
      do while (feh.lt.zmin .or. feh.gt.zmax)
         X = ran3(seed)
         i = floor(X*nzdist) + 1
         feh = zdist(i)
         feh = feh + (ran3(seed) - 0.5d0) * zfuzz
      end do

c convert from [Fe/H] to feh
      getzcum = 2.d-2 * 10.d0**feh

      return
      end



