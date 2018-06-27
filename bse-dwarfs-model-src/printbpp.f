      subroutine printbpp
      implicit none
      INTEGER j, k, kw, kstar(2)
      REAL bcm(50000,34),bpp(80,10)
      COMMON /BINARY/ bcm,bpp

      CHARACTER*8 label(14)

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
      if(bpp(j,10).le.0.0) goto 60
      kstar(1) = INT(bpp(j,4))
      kstar(2) = INT(bpp(j,5))
      kw = INT(bpp(j,10))
      WRITE(*,99100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),label(kw)
      goto 52
 60   continue
99100 FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)
      WRITE(*,*)
 100  FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)
      RETURN
      END
