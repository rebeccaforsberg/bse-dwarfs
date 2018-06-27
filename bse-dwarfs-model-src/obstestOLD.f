      logical function obstest2OLD(mbol1,mags1,mbol2,mags2,Dgal,sep,istar)
      implicit none
      real*8 mbol1,mags1(*),Dgal,mbol2,mags2(*),sep
      integer istar

      include "rand.h"
      include "ubvgrid.h"
      include "obscount.h"

      real*8 absMg1,absmg2,mu,Mg,mags(nbands),mg1,mg2,angsep,frat,mbol,f
      integer i

      real*8 aursun
      parameter (aursun=214.95d0)

      obstest2OLD = .false.

c Absolute G-band magnitudes
      absMG1 = mbol1 - mags1(2)
      absMG2 = mbol2 - mags2(2)

c Distance modulus (Dgal in kpc)
      mu = 5.d0 * log10(Dgal/0.01d0)

c Apparent G-band magnitudes
      MG1 = absMG1 + mu
      MG2 = absMG2 + mu

c G-band flux ratio
      frat = 1.d1**(-0.4d0*(MG1-MG2))

c Angular separation in arcsec
      angsep = sep/aursun/(1.d3*Dgal)

c Large separation: both observed
      if (angsep.gt.2.d0) then
         if (MG1.le.11.5d0) then
            call writeSingleStar(mbol1, mags1, dgal, mu, istar)
            nobs = nobs + 1
            obstest2old = .true.
         end if
         if (MG2.le.11.5d0) then
            call writeSingleStar(mbol2, mags2, dgal, mu, istar)
            nobs = nobs + 1
            obstest2old = .true.
         end if
         return
      end if

c Blend light for small angular separations or large flux ratios
      if (angsep.gt.0.1d0) then
         if (frat.lt.1.d1 .and. frat.gt.0.1d0) then
            return
         end if
      end if

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
      MG = mbol - mags(2) + mu

c Complete to 11.5
      if (MG .gt. 11.5d0) return
      
      call writeSingleStar(mbol, mags, dgal, mu, istar)
      nobs = nobs + 1
      obstest2old = .true.

      return
      end


c Write out data for an observed single star
      subroutine writeSingleStar(mbol, mags, Dgal, mu, istar)
      implicit none
      real*8 mbol,mags(*),Dgal,mu
      integer istar

c Defined fixed magnitude error
      real*8 magerr
      parameter (magerr=0.025d0)

      real*8 obs(8,2),absMg,Mg,plx,plxerr,f
      integer jbin,i
      real*8 gasdev

      integer nplxbins
      parameter (nplxbins = 7)
      real*8 plxmags(nplxbins), plxerrs(nplxbins)
      data plxmags, plxerrs /6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
     :                       399., 348., 327., 407., 601., 722., 702./

      include "rand.h"
      include "ubvgrid.h"

c Absolute Gaia G-band magnitude
      absMg = mbol - mags(4)

c Apparent G-band magnitude
      Mg = absMg + mu

c Parallax in mas
      plx = 1.d0 / Dgal

c Parallax error
      if (Mg.lt.plxmags(1)) then
         plxerr = plxerrs(1)
      else if (Mg.gt.plxmags(nplxbins)) then
         plxerr = plxerrs(nplxbins)
      else
c The next line RELIES on the difference between mag bins being 1 each
         jbin = floor(Mg-plxmags(1)) + 1
         f = Mg-plxmags(jbin)
         if (f.gt.1.d0 .or. f.lt.0.d0) then
            write (*,*) 'FUCKUP: ', jbin, plxmags(jbin), plxmags(jbin+1), Mg
            stop
         end if
         plxerr = f * plxerrs(jbin+1) + (1.d0 - f) * plxerrs(jbin)
      end if

c Apply the "observational" uncertainties
      obs(1,1) = plx 
      obs(1,2) = plx + plxerr * gasdev() / 1.d3
      obs(2,1) = Mg
      obs(2,2) = Mg  + magerr * gasdev()
      do i = 1, 6
         obs(i+2,1) = mbol - mags(i) 
         obs(i+2,2) = obs(i+2,1) + magerr * gasdev()
      end do

      write (30, 99001) (obs(i,1),i=1,8),istar
      write (29, 99001) (obs(i,2),i=1,8),istar
99001 FORMAT(1p, 8(1x,e10.3),i9)

      return
      end








