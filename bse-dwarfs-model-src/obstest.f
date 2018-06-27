c     subroutine obstest2(mbol1,mags1,mbol2,mags2,Dgal,sep,istar,ota,otb,obt2)
      subroutine obstest2(mbol1,mags1,logg1,logt1,mbol2,mags2,logg2,logt2,Dgal,
     :                           sep,istar,ota,otb,obt2,absMJ,absMK,mbol,MH,frat,
     :                           kw,kw2)
      implicit none
      real*8 mbol1,mags1(*),Dgal,mbol2,mags2(*),sep,logg1,logt1,logg2,logt2
      integer istar,kw,kw2
      logical ota,otb,obt2,obstest1

      include "rand.h"
      include "ubvgrid.h"
      include "obscount.h"

      real*8 absMH1,absmH2,mu,MH,mags(nbands),mH1,mH2,angsep,frat,mbol,f,absMJ,absMK,JKS
      integer i

      real*8 aursun
      parameter (aursun=214.95d0)

      obt2 = .false.
      ota = .false.
      otb = .false.

c Absolute H-band magnitudes
      absMH1 = mbol1 - mags1(2)
      absMH2 = mbol2 - mags2(2)

c Distance modulus (Dgal in kpc)
      mu = 5.d0 * log10(Dgal/0.01d0)

c Apparent H-band magnitudes
      MH1 = absMH1 + mu
      MH2 = absMH2 + mu

c "Apparent" bolometric magnitude

      MH1 = mbol1 + mu
      MH2 = mbol2 + mu


c Combine fluxes
      do i = 1, nbands
         f = 0.d0
         if (mags1(i).gt.-4.d1) f = f + 1.d1**(-.4d0*mags1(i))
         if (mags2(i).gt.-4.d1) f = f + 1.d1**(-.4d0*mags2(i)) 
         if (f.gt.0.d0) then 
            mags(i) = -2.5d0*log10(f)
         else
            mags(i) = -4.d1
         end if
      end do
      f = 0.d0
      if (mbol1.gt.-4.d1) f = f + 1.d1**(-.4d0*mbol1)
      if (mbol2.gt.-4.d1) f = f + 1.d1**(-.4d0*mbol2) 
      mbol =  -2.5d0*log10(f)
      MH = mbol - mags(2) + mu

c Make "apparent" bolometric magnitude 

      MH = mbol + mu

c If the combined binary is too faint to see then we are wasting our time
      if (MH.gt.20) return

c Angular separation in arcsec (sep should be projected onto sky)
      angsep = sep/aursun/(1.d3*Dgal)

c Large separation: both potentially observed
      if (angsep.gt.6.d0) then
c Do tests on both stars individually to see whether they are bright enough to be seen
         ota = obstest1(mbol1,mags1,Dgal,logg1,istar,absMJ,absMK,MH)
         otb = obstest1(mbol2,mags2,Dgal,logg2,istar+1,absMJ,absMK,MH)
         return
      end if

c Colour check (if binary is too red we will see nothing)
      absMJ = mbol - mags(1)
      absMK = mbol - mags(3)
      JKs = absMJ - absMK
c      if (JKs.lt.0.5d0) return

c H-band flux ratio to work out which star is brighter
      frat = 1.d1**(-0.4d0*(MH1-MH2))

      if (frat.gt.1.d1) then
c Only see star 1
         otb = .false.
         ota = .true.
c Check for star 1 being a giant
         if (logg1.le.4.0d0) then
            ota = .false.
         else
            ota = .true.
         end if
      else if (frat.lt.1.d-1) then
c Only see star 2
         ota = .false.
         otb = .true.
         if (logg2.le.4.0d0) then
            otb = .false.
         else
            otb = .true.
         end if
      else
c SB2: check that they aren't both giants!
      ota = .true.
      otb = .true.
         if (logg1.gt.4.0d0 .and. logg2.gt.4.0d0) then
            ota = .true.
            otb = .true.
         else
            ota = .false.
            otb = .false.
         end if
      end if

c Do we observe the binary, somehow?  This is .false. for observing separate stars from the binary
      if (ota.or.otb) then
         obt2 = .true.
         nobs = nobs + 2
         seen = .true.
c      else if ((kw.gt.1 .and. kw.lt.10) .or. (kw2.gt.1 .and. kw2.lt.10)) then
cc XXX TEMP keep all the systems where we have at least one plausible star REGARDLESS
c         obt2 = .true.
c         ota = .true.
c         otb = .true.
c         nobs = nobs + 2
c         seen = .true.
      end if

      return
      end

c Test whether a single star will be observed and if so observe it
      logical function obstest1(mbol,mags,Dgal,logg,istar,absMJ,absMK,MH)
      implicit none
      real*8 mbol,mags(*),Dgal,logg
      integer istar

      include "rand.h"
      include "ubvgrid.h"
      include "obscount.h"

      real*8 absMH,mu,MH,JKs,absMJ,absMK

      obstest1 = .false.

c Absolute H-band magnitude
      absMJ = mbol - mags(1)
      absMH = mbol - mags(2)
      absMK = mbol - mags(3)

c Distance modulus (Dgal in kpc)
      mu = 5.d0 * log10(Dgal/0.01d0)

c Apparent H-band magnitude
      MH = absMH + mu

c Redefine as "apaprent" bolometric magnitude

      MH = mbol + mu      

c Complete to ~20
      if (MH .gt. 20d0) return

c Cut on J-K_s colour
      JKs = absMJ - absMK
c      if (JKs.lt.0.5d0) return

c Minimum log g to avoid giants
      if (logg.le.4.0d0) return
      
      nobs = nobs + 1
      seen = .true.
      obstest1 = .true.

      return
      end
