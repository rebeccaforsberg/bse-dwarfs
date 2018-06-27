***
      PROGRAM bse
***
*
* Evolves a binary by calling evolv2.f 
* (see header of subroutine for algorithm description). 
*
* Required input is described below. 
***
* See Tout et al., MNRAS, 1997, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code.
* Updated reference is:
*           Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
* (please use this one).
***
* For single star evolution see Hurley, Pols & Tout, 2000, MNRAS, 315, 543.
* or Hurley, 2000, PhD Thesis, University of Cambridge (Chapter 2).
* The binary evolution algorithm is described in Chapter 3 of the thesis.
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
***
      implicit none
*
      INCLUDE 'const_bse.h'
      INCLUDE 'rand.h'
      INCLUDE 'ubvgrid.h'
      INCLUDE 'obscount.h'
*
      integer nstarmax,apogee_starmax,japogee,napogee
      parameter(nstarmax=10000000)
      parameter(apogee_starmax=1000000)
      integer kw,kw2,kstar(2),j,k,time,i,n,ofile,status1,status2,ntimes,l,m,istar
      integer nignored,nunseens,nunseenb

c Input
      integer inum,nsys,imfm1,imfm2,ibin,imet,iage,igal,noff,iwhich
      real*8  mmin1,mmax1,mmin2,mmax2,gammaq,fbin,zmin,zmax,zfuzz,agemin,agemax,tauage

      real*8  xmin,xmax,Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal,m1,m2,X
      real*8  getmass,dm91,r10,getzcum

*   
      real*8 logl1,logt1,mass1,logz1,logg1,logm1,t1,logr1,t
      real*8 bolc1,mv1,uminb1,bminv1,vmini1,logg2,logm2,t2,logr2
      real*8 logl2,logt2, mass2,logz2,magsun1,mags1(nbands)
      real*8 bolc2,mv2,uminb2,bminv2,vmini2,magsun2,mags2(nbands)
      real*8 mass0(2),mass(2),z,zpars(20),logmh
      real*8 epoch(2),tms(2),tphys,tphysf,dtp,aj
      real*8 rad(2),lum(2),ospin(2),times(28)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 tb,ecc,yearsc,sep,xdiffs(28),vdiffs(28),meanv,sdv,fnt
      real*8 indata(6,nstarmax),gconst,indata2(29,apogee_starmax)
      PARAMETER(yearsc=3.1557d+07)
      parameter(gconst=-10.6071d0)
      CHARACTER*8 label(14)
      CHARACTER*140 line
      CHARACTER*50 fn
      real*8 ran3
      logical single,unfinished,obstest1,ot2,ota,otb,galReject

*
************************************************************************
* Input:
*
* mass is in solar units.
* tphysf is the maximum evolution time in Myr.
* tb is the orbital period in days.
* kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state. 
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* eccentricity can be anywhere in the range 0.0 -> 1.0.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). #defunct#
* ceflag = 3 activates de Kool common-envelope model (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*
* If you enter a negative kstar then parameters for an evolved star are
* required in the order of:
* current age, initial mass and spin rate, 
* otherwise the star will start on the ZAMS.
*

c Tables of data
      CALL loadtab

      kstar(1) = 1
      kstar(2) = 1
      OPEN(22,file='binary.in', status='old')
      READ(22,*)neta,bwind,hewind,alpha1,lambda
      READ(22,*)ceflag,tflag,ifflag,wdflag,bhflag,nsflag,mxns,idum
      READ(22,*)pts1,pts2,pts3
      READ(22,*)sigma,beta,xi,acc2,epsnov,eddfac,gamma
      if(kstar(1).lt.0.or.kstar(2).lt.0)then
         READ(22,*)tphys
         READ(22,*)aj,mass(1),ospin(1)
         epoch(1) = tphys - aj
         kstar(1) = ABS(kstar(1))
         READ(22,*)aj,mass(2),ospin(2)
         epoch(2) = tphys - aj
         kstar(2) = ABS(kstar(2))
      else
*
* Initialize the parameters.
* Set the initial spin of the stars. If ospin is zero (actually < 0.001)
* at time zero then evolv2 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
* If ospin is negative then the stars will be in co-rotation with the orbit.
*
         tphys = 0.d0
         mass(1) = mass0(1)
         epoch(1) = 0.d0
         ospin(1) = 0.d0
         mass(2) = mass0(2)
cc         mass(2) = 0.d0
         epoch(2) = 0.d0
         ospin(2) = 0.d0
      endif
      if(idum.gt.0) idum = -idum
      CLOSE(21)
      CLOSE(22)
      WRITE(*,*)
*
* Note that this routine can be used to evolve a single star if you 
* simply set mass(2) = 0.0 or tb = 0.0 (setting both is advised as  
* well as some dummy value for ecc). 
*
************************************************************************
*
* Set parameters which depend on the metallicity 
*      write (*,*) 'Going into zcnsts, z=', z

*      CALL zcnsts(z,zpars)
      
*      write (*,*) 'Back from zcnsts'
*
* Set the collision matrix.
*
      
      CALL instar
*
      label(1) = 'INITIAL '
      label(2) = 'KW CHNGE'
      label(3) = 'BEG RCHE'
      label(4) = 'END RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO REMNT'
      label(10) = 'MAX TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG SYMB'
      label(13) = 'END SYMB'
      label(14) = 'BEG BSS'
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the bcm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = 0.d0
      z = 0.02d0
*
* Output files
      unfinished = .true.
      m1 = 0.8d0
      istar = 1
      do while (unfinished)
         mass0(1) = m1
         mass(1) = mass0(1)

         write (fn,99900) m1
99900 FORMAT("track.", f6.3, ".out")
         open (26, file=fn, status="new")

c Deal with single star case
         single = .true.
         mass0(2) = 0.d0
         mass(2) = 0.d0
         tb = 0.d0
         ecc = 0.d0

         call zcnsts(z,zpars)

c Other BSE setup
         kstar(1) = 1
         kstar(2) = 1
         tphys = 0.d0
         mass(1) = mass0(1)
         epoch(1) = 0.d0
         ospin(1) = 0.d0
         mass(2) = mass0(2)
cc         mass(2) = 0.d0
         epoch(2) = 0.d0
         ospin(2) = 0.d0

         tphysf = 1.4d4    ! In Myr

c Write out initial stellar properties
c         write (line,99001) istar, mass0(1), mass0(2), tb, ecc, z, tphysf, 
c     &                    Rgal,qgal,zgal,Dgal,lgal,bgal,xgal,ygal
         
         CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &              menv,renv,ospin,epoch,tms,
     &              tphys,tphysf,dtp,z,zpars,tb,ecc)

c Check for some results and write output
         j = 1
         write (*,*) (bcm(j,i),i=1,20)
         do while (bcm(j,1).ge.0.d0)
            
            kw = INT(bcm(j,2))
            kw2 = INT(bcm(j,16))
           
c Assemble stellar parameters
            mass1 = bcm(j,4)
            if (mass1 .gt. 0.d0) then
               logm1 = log10(mass1)
            else
               logm1 = -50.d0
            end if
            logl1 = bcm(j,5)
            logr1 = bcm(j,6)
            logt1 = bcm(j,7)
            logg1 = logm1 - 2.d0*logr1 + 4.4383174
            magsun1 = 4.75d0 - 2.5d0*logl1                  

            mass2 = bcm(j,18)
            if (mass2 .gt. 0.d0) then
               logm2 = log10(mass2)
            else
               logm2 = -50.d0
            end if
            logl2 = bcm(j,19)
            logr2 = bcm(j,20)
            logt2 = bcm(j,21)
            logg2 = logm2 - 2.d0*logr2 + 4.4383174                  
            magsun2 = 4.75d0 - 2.5d0*logl2                  

c Metallicity
            logmh = log10(z/0.02d0)

c Initially single star
            CALL gt2mag(logg1,logt1,logmh,mags1,status1)
            write (26,999) istar,bcm(j,1),kw,kw2,bcm(j,4),
     &                bcm(j,8),bcm(j,6),
     &                bcm(j,5),
     &                bcm(j,7),bcm(j,13),
     &                logg1,
     &                magsun1,(magsun1-mags1(n),n=1,nbands),
     &                status1   
            j = j + 1
         end do
         close(26)
         m1 = m1 + 0.10d0
         if (m1.gt.1.d1) unfinished = .false.
         istar = istar + 1
      end do

 80   FORMAT('#istar    TIME       kw kw2 mass1    mass2    massc1   massc2   ',    ! 8
     &       'logr1      logr2      radrol1    radrol2    ',                        ! 12
     &       'logl1      logl2      logt1      logt2      ',                        ! 16
     &       'ospin1      ospin2      sep         tb          ecc         ',        ! 21
     &       'logg1           logg2           magsun1         magsun2         ',    ! 25
     &       'J1              H1              Ks1             J2              ',    ! 29
     &       'H2              Ks2            ',                                     ! 31
     &       'stat1 stat2 ntims ',                                                  ! 34
     &       'meanv       stddevv     ',                                            ! 36
     &       'xdiffs1     xdiffs2     xdiffs3     xdiffs4     ',                    ! 40
     &       'xdiffs5     xdiffs6     xdiffs7     xdiffs8     ',
     &       'xdiffs9     xdiffs10    xdiffs11    xdiffs12    ',
     &       'xdiffs13    xdiffs14    xdiffs15    xdiffs16    ',
     &       'xdiffs17    xdiffs18    xdiffs19    xdiffs20    ',
     &       'xdiffs21    xdiffs22    xdiffs23    xdiffs24    ',
     &       'xdiffs25    xdiffs26    xdiffs27    xdiffs28    ',
     &       'vdiffs1     vdiffs2     vdiffs3     vdiffs4     ',
     &       'vdiffs5     vdiffs6     vdiffs7     vdiffs8     ',
     &       'vdiffs9     vdiffs10    vdiffs11    vdiffs12    ',
     &       'vdiffs13    vdiffs14    vdiffs15    vdiffs16    ',
     &       'vdiffs17    vdiffs18    vdiffs19    vdiffs20    ',
     &       'vdiffs21    vdiffs22    vdiffs23    vdiffs24    ',
     &       'vdiffs25    vdiffs26    vdiffs27    vdiffs28    ',
     &       'iw')                                                                  ! 93


 82   FORMAT('#TIME kw kw2 mass1   mass2   massc1   massc2',
     &       '   logr1   logr2   radrol1   radrol2   logl1   logl2',
     &       '   logt1   logt2   ospin1   ospin2   sep   tb   ecc',  
     &       '   logg1 logg2 magsun1 magsun2 J1 H1 Ks1',
     &       '  J2 H2 Ks2 status1 status2')     
     
 81   FORMAT('#TIME  kw  kw2  mass  massc',
     &       '   logr  logl',
     &       '   logt  ospin  logg  magsun  J  H  Ks status')
c     &       '  I  J  H  K  status')

   83 FORMAT('#istar        m1       m2     period     ecc    z      age        ',
     :       '  Rgal     qgal     zgal     Dgal     lgal     bgal     xgal     ygal')

 9    FORMAT(i9,1x,f10.4,2i3,10(5x,f10.4),7(5x,e15.4),10f15.4,2i6)
 99   FORMAT(i9,1x,(f10.4),2(i3),4(1x,f8.4),8(1x,f10.4),1P,5(1x,e11.4),0P,10(1x,f15.4),
     &       3(i6),1P,58(1x,E11.4),0P,1x,i2)
 999  FORMAT(i9,1x,f11.4,2i3,8(1x,f10.4),3(1x,f10.4),1(1x,i4))
 
99001 FORMAT(i9,1x,2(f8.3,1x),1p,e10.3,0p,1x,f6.3,1x,f6.3,1x,1p,e10.3,0p,
     :       8(1x,f8.3))
99002 FORMAT('#istar    mass0(1) mass0(2) tb         ecc    z      ',
     :       'tphysf     ',
     :       'Rgal     qgal     zgal     Dgal     lgal     bgal     xgal     ygal     ')
      
*
************************************************************************
* Output:
* First check that bcm is not empty.
*
*      if(bcm(1,1).lt.0.0) goto 50
*
* The bcm array stores the stellar and orbital parameters at the 
* specified output times. The parameters are (in order of storage):
*
*    Time, 
*    [stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, core radius, mass of any convective 
*    envelope, radius of the envelope, epoch, spin, mass loss rate and 
*    ratio of radius to roche lobe radius (repeated for secondary)],
*    period, separation, eccentricity.
*
 
*     OPEN(23,file='binary.dat', status='unknown')
*      j = 0
* 30   j = j + 1
*      if(bcm(j,1).lt.0.0)then
*         bcm(j-1,1) = bcm(j,1)
*         j = j - 1
*      endif
*      kw = INT(bcm(j,2))
*      kw2 = INT(bcm(j,16))
*      WRITE(23,99)bcm(j,1),kw,kw2,bcm(j,4),bcm(j,18),
*     &            bcm(j,8),bcm(j,22), 
*     &            bcm(j,6),bcm(j,20),bcm(j,15),bcm(j,29),
*     &            bcm(j,5),bcm(j,19),bcm(j,13),bcm(j,27),
*     &            bcm(j,14),bcm(j,28),
*     &            bcm(j,31),bcm(j,32)
*      if(bcm(j,1).ge.0.0) goto 30
*      CLOSE(23)
*
* The bpp array acts as a log, storing parameters at each change
* of evolution stage.
*
 50   j = 0
      WRITE(*,*)'     TIME      M1       M2   K1 K2        SEP    ECC',  
     &          '  R1/ROL1 R2/ROL2  TYPE'
 52   j = j + 1
      if(bpp(j,1).lt.0.0) goto 60
      kstar(1) = INT(bpp(j,4))
      kstar(2) = INT(bpp(j,5))
      kw = INT(bpp(j,10))
      WRITE(*,99100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),label(kw)
      goto 52
 60   continue
99100 FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)
      WRITE(*,*)
*
************************************************************************
*
      STOP
      END
***
