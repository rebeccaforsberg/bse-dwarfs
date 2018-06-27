      SUBROUTINE loadTab
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"

      CHARACTER*2048 line
      CHARACTER*256 filenames(Nmh)

      DATA filenames / "bctab_m25.dat", "bctab_m20.dat", "bctab_m15.dat",  
     :                 "bctab_m10.dat", "bctab_m05.dat", "bctab_p00.dat",  
     :                 "bctab_p05.dat"/
      DATA mhs       / -2.5, -2, -1.5, -1, -.5, 0., 0.5 /

      INTEGER imh,i
      LOGICAL lgotgrid

      call blankGrid

      lGotGrid = .false.
c Loop over the metalicities
      do imh = 1, nmh
         write (*,*) 'Loading metallicity ', mhs(imh), ' from file ', filenames(imh)
         open(20,file=filenames(imh),status="old",form="formatted")
         call loadGridFromOpenFile(imh,lgotgrid)
         close(20)
      end do

c Giant grid is back to front
      call reverseGiantGrid

c Convert the grid in Teff to log10(teff)
      do i = 1, nlogt
         logts(i) = log10(logts(i))
      end do

c Obtain the minimum temeperature for which we have normal giant models
      jTminGiant = 1
      do while (grid(1,jTminGiant,1,1).lt.-1e19)
         jTminGiant = jTminGiant + 1
      end do
      logTminGiant = logts(jTminGiant)
      return
      end

c Reverse the order of entries in the M-giant grid (sigh)
      subroutine reverseGiantGrid
      implicit none
      include "ubvgrid.h"
      
      integer i,j,k,imh
      real*8 tmp(nbands)

      do i = 1, ngiant/2
         j = ngiant+1-i
         if (i<j) then
            tmp(1) = giantts(i)
            giantts(i) = giantts(j)
            giantts(j) = tmp(1)
            do imh = 1, nmh
               do k = 1, nbands
                  tmp(k) = giantGrid(k,i,imh)
                  giantGrid(k,i,imh) = giantGrid(k,j,imh)
                  giantGrid(k,j,imh) = tmp(k)
               end do
            end do
         end if
      end do
      return
      end
            

c This routine does the work.  Note that teff is unlogged in the input blob.
      SUBROUTINE loadGridFromOpenFile(imh,lgotgrid)
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"

      INTEGER imh
      LOGICAL lgotgrid

      CHARACTER*2048 line
      INTEGER nl,mode,i,il,ix,iy,j
      REAL*8 teff,logg,bc(nbands)
      REAL*8 inlines(nbands+2,nlinesmax)
      INTEGER bisect

      open (13,file="tableDebug.out",access="append")

      nl = 0
      mode = 0
      lgotgrid = .false.
c The following two counters are re-set for each file
      nhot = 0
      ngiant = 0
      do
         read (20,'(A2048)',end=100) line
         if (line(1:1).eq."#") then
c If we have just started or we are already in the M-giant bit ignore comments
            if (nl.eq.0.or.mode.eq.2) then
               cycle
            else
c Else a comment suggests the start of the M-giant bit 
               mode = 2
               cycle
            end if
         end if

         read (line,*) il, teff, logg, (bc(i),i=1,nbands)

         if (mode.eq.0) then
cc If we have an odd temperature switch to the hot grid
c            if (mod(teff,10.d0).ne.0) then
c               mode = 1
c            else
               nl = nl + 1
               inlines(1,nl) = teff
               inlines(2,nl) = logg
               do i = 1, nbands
                  inlines(i+2,nl) = bc(i)
               end do
               cycle
c            end if
         end if

cc Mode==1 => hot grid
c         if (mode.eq.1) then
c            nhot = nhot + 1
c            hotTs(nhot) = log10(teff)
c            do i = 1, nbands
c               hotgrid(i,nhot,imh) = bc(i)
c            end do
c            cycle

c Mode==2 => M-giant bit
         if (mode.eq.2) then
            ngiant = ngiant + 1
            giantTs(ngiant) = log10(teff)
            do i = 1, nbands
               giantgrid(i,ngiant,imh) = bc(i)
            end do
            cycle
         end if

         write (*,*) 'ERROR: CONFUSED in loadTab', mode
         stop
      end do

c We get to here once we reach the end of the file
100   continue
      write (*,*) 'Giant grid: ngiant=', ngiant
      write (*,*) (giantgrid(1,i,imh),i=1,min(ngiant,10))
      write (*,*) '...'

c Now sort out the grid
      if (.not.lgotgrid) then
         call extractGridEdge(inlines,logts,1,nl,nlogt)
         call extractGridEdge(inlines,loggs,2,nl,nlogg)
         lgotgrid = .true.
         write (*,*) 'Extracted logg and logt values'
         do i = 1, min(10,nlogt,nlogg)
            write (*,*) i, logts(i), loggs(i)
         end do
         write (*,*) '...'
      end if

c Stuff the input data blob into the grid
      do i = 1, nl
         teff = inlines(1,i)
         logg = inlines(2,i)
         ix = bisect(teff,logts,nlogt)
         iy = bisect(logg,loggs,nlogg)
         do j = 1, nbands
            grid(j,ix,iy,imh) = inlines(j+2,i)
         end do
      end do
      write (*,*) 'Inserted blob into grid', imh
      do i = 1, min(10,nlogt,nlogg)
         write (*,*) i, logts(i), loggs(i)
      end do

c Infill the gaps for which data is not provided
      call infill(imh)

      close(13)

      return
      end

c Fill in "holes" in the grid by interpolating along the temperature axis
      SUBROUTINE infill(imh)
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"
      INTEGER imh

      INTEGER ix,iy,ilast,i,j
      REAL*8 ax(nlogtmax),x,omx,x1,x2
      logical found

c Cache a set of logged logt values (YA RLY)
      do iy = 1, nlogt
         ax(iy) = dlog10(logts(iy))
      end do

c Loop over log g values
      do ix = 1, nlogg

c Find the first log t value for this logg value
         ilast = 0
         do iy = 1, nlogt
            if (grid(1,iy,ix,imh).gt.-1.d19) then
               ilast = iy
               write (13,*) imh, ix, iy, loggs(ix), logts(iy), 1
               exit
            else
               write (13,*) imh, ix, iy, loggs(ix), logts(iy), 0
            end if
         end do
         if (ilast.eq.0) then
            write (*,*) 'ERROR: FAILED to find any values in infill', 
     :                       ix, loggs(ix), imh, mhs(imh)
            stop
         end if
         mint(ix) = ax(ilast)
         maxt(ix) = ax(ilast)

c Now loop through all the others
         do iy = ilast+1, nlogt
            if (grid(1,iy,ix,imh).gt.-1.d19) then
               maxt(ix) = ax(iy)
c There is a value here
               if (iy-ilast.eq.1) then
c There is no gap; shuffle ilast up and carry on
                  ilast = iy
                  write (13,*) imh, ix, iy, loggs(ix), logts(iy), 1
               else
c We have found both ends of a gap; infill
                  x1 = ax(ilast)
                  x2 = ax(iy)
                  write (13,*) imh, ix, iy, loggs(ix), logts(iy), 1
                  do i = ilast+1, iy-1
                     x = (ax(i)-x1) / (x2-x1)
                     omx = 1.d0 - x
                     do j = 1, nbands
                        grid(j,i,ix,imh) = omx*grid(j,ilast,ix,imh) + x*grid(j,iy,ix,imh)
                     end do
                     write (13,*) imh, ix, i, loggs(ix), logts(i), 2
                  end do
                  ilast = iy
               end if
            else
c There is not a value here
               i=iy   ! NOP
               write (13,*) imh, ix, i, loggs(ix), logts(i), 3
            end if
         end do
      end do

c Loop over all logts to find the min and max logg for each point
      do ix = 1, nlogt
         found = .false.
         do iy = 1, nlogg
            if (grid(1,ix,iy,imh).gt.-1.d19) then
               if (found) then
                  maxg(ix) = loggs(iy)
               else
                  ming(ix) = loggs(iy)
                  maxg(ix) = loggs(iy)
                  found = .true.
               end if
            end if
         end do
      end do
            
      return
      end
        

c Find the i in [1:nx] that minimises |x(i)-y| where x is a sorted array
      INTEGER FUNCTION bisect(y,x,nx)
      IMPLICIT NONE
      real*8 y, x(*)
      integer nx

      integer i1,i2,i
      real*8 d1, d2

      i1 = 1
      i2 = nx
      do while (i2-i1.gt.1)
         i = (i1+i2)/2
         if (x(i).lt.y) then
            i1 = i
         else
            i2 = i
         end if
      end do
      d1 = dabs(y-x(i1))
      d2 = dabs(y-x(i2))
      if (d1.lt.d2) then
         bisect = i1
      else
         bisect = i2
      end if
      return
      end

            
c Extract the grid edge from the input data blob
      SUBROUTINE extractGridEdge(inlines,outData,ix,nl,nOutData)
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"
      REAL*8 inlines(nbands+2,nlinesmax), outdata(*)
      REAL*8 dum(nlinesmax)
      INTEGER ix, nl, noutdata

      INTEGER i,j

      write (*,*) 'Extracting data with ix=', ix, ' and nl=', nl
      
c Take a copy of the data to be sorted
      do i = 1, nl
         dum(i) = inlines(ix,i)
      end do

c Sort it
      call sort(nl,dum)

c Uniq it
      nOutData = 1
      outData(1) = dum(1)
      do i = 2, nl
         if (dum(i).ne.outData(nOutData)) then
            nOutData = nOutData + 1
            outData(nOutData) = dum(i)
         end if
      end do

      write (*,*) 'Extracted ', noutdata, ' lines'
      do i = 1, min(noutdata,10)
         write(*,*) i, outdata(i)
      end do
      write (*,*) '...'

      return
      end




      






      SUBROUTINE blankGrid
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"

      INTEGER i, j, k, l

c Set the main grid to the default (unfilled) value
      DO i = 1, nmh
         do j = 1, nloggmax
            do k = 1, nlogtmax
               do l = 1, nbands
                  grid(l,k,j,i) = -1.d20
               end do
            end do
         end do
      end do
      nhot = 0
      ngiant = 0

      return
      end





