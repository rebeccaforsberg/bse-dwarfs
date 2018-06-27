   

      SUBROUTINE countEvents(eventCounts, rlofEventString)
      IMPLICIT NONE
      INCLUDE "const_bse.h"
      INTEGER eventCounts(3)
      CHARACTER*30 rlofEventString

      INTEGER i,iev,Nrlof,cee,kw1,kw2,starrlofing
      INTEGER rloflist(6)
      real*8 hestarmass
      LOGICAL inrlof


      inrlof = .FALSE.
      Nrlof = 0
      cee = 0;
      hestarmass = 0.d0

c Blank counters of RLOF events
      do i = 1, 3
         eventCounts(i) = 0
      end do
      do i = 1, 6
         rloflist(i) = 0
      end do
      do i = 1, 80
         kw1 = int(bpp(i,4))
         kw2 = int(bpp(i,5))
         iev = int(bpp(i,10))
c Break out of the loop if the next datum is an endpoint / incomprehensible
         if (bpp(i,1).lt.0.d0 .or. iev.le.0.0) exit

c Check for coelescence
         if (iev.eq.6) eventcounts(3) = eventcounts(3) + 1
   
c         /* Check for RLOF start/end */
         if (bpp(i,10).eq.3.or.(bpp(i,10).eq.5.and.(.not.inrlof))) then
            if (inrlof) then
               write (*,*) "WARNING!!!!!!!! inrlof = ", inrlof
            end if
            inrlof = .true.
            Nrlof = Nrlof + 1
            if (Nrlof.ge.5) then
               write(*,*) "Lots of RLOF filling events!"
               Nrlof = 5
            end if
c            /* Which star is RLOFing? 0 == primary, 1 == secondary 
            if (bpp(i,8).gt.bpp(i,9)) then
               starRlofing = 0 
            else 
               starRlofing = 1
            end if
            rloflist(Nrlof) =rlofList(Nrlof) + starRlofing
c            /* Type of two stars */
            rloflist(Nrlof) = rloflist(Nrlof) + 
     :                        int(bpp(i,4+starRlofing)) * 2
            rloflist(Nrlof) = rloflist(Nrlof) + 
     :                        int(bpp(i,4+1-starRlofing)) * 32
            cee = 0
            cycle
         end if
c         /* End of RLOF */
         if (bpp(i,10) .eq. 4) then
            inrlof = .false.
            if (nrlof.le.0) then
               write (*,*) 'NRLOF=', NRLOF
               eventcounts(2) = -1
               return
            end if
            rloflist(Nrlof) = rloflist(Nrlof) + cee * 512
            cycle
         end if
c         /* Check for CEE */
         if (bpp(i,10) .eq. 7) then
            eventcounts(2) = eventcounts(2) + 1
            cee = 1
c            /* Check for CEE following CONTACT */
            if (bpp(i-1,10) .eq. 5) then
               inrlof = .false.
               rloflist(Nrlof) = rloflist(Nrlof) + cee * 512
            end if
            cycle
         end if
      end do

      eventcounts(1) = nrlof
   
c       Create a nice string 
      write(rlofeventstring,99001) (min(9999,rloflist(i)),i=1,6)
99001 format(6(i4,1x),f6.3) 
      return
      end

