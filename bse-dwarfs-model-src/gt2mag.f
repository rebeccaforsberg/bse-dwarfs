      SUBROUTINE gt2mag(logg,logt,logmh,mags,status)
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"

      REAL*8 logg,logt,logmh,mags(*)
      INTEGER status

      real*8 zlm,zlt,zlg
      integer tag(2,2),i,j,k,l,ix,iy,iz,ix1,ix2,iy1,iy2,iz1,iz2,ttot,off(2),kkk,
     :        usegiant,testinmgiant
      real*8 f(3),fx,fy,fz,d2,magtmp(Nbands),fmix
      common /spy/ ix1,iy1,iz1

      logical debug

      status = 0
      usegiant = 0
      fmix = 0.d0

      do i = 1, nbands
         mags(i) = 0.d0
      end do

c First check to see whether we are off the metallicity scale
      if (logmh.lt.mhs(1)) then
         status = 1
         zlm = mhs(1)
      else if (logmh.gt.mhs(nmh)) then
         status = 1
         zlm = mhs(nmh)
      else
         zlm = logmh
      end if
      call bisect2(zlm,mhs,nmh,iz1,iz2)
      fz = (zlm-mhs(iz1))/(mhs(iz2)-mhs(iz1))

c Special check for giants, which need to be treated carefully
      if (logg.lt.1.75d0) then
         call selectInterpolationGiant(logg, logt,
     :         zlg, zlt, ix1, ix2, iy1, iy2, debug)
         status = status + 2048
         goto 100
      end if

c Check for being off the logg and logt scales
      if (logt.gt.logts(nlogt)) then
         off(1) = 1
      else if (logt.lt.logts(1)) then
         off(1) = -1
      else
         off(1) = 0
         zlt = logt
      end if
      if (logg.gt.loggs(nlogg)) then
         off(2) = 1
      else if (logg.lt.loggs(1)) then
         off(2) = -1
      else
         off(2) = 0
         zlg = logg
      end if

c Off both scales: nasty!
      if (off(1).ne.0.and.off(2).ne.0) then
         call findNearestPoint(logt, logg, fz, iz1, iz2, mags)
         status = status + 2
         return
      end if

c Off top of log g scale (A): use nearest log g and interpolate in log t
      if (off(2).eq.1) then
         call bisect2(zlt,logts,nlogt,ix1,ix2)
         zlg = min(maxg(ix1),maxg(ix2))
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
         status = status + 4
         go to 100
      end if

c Off bottom of log g scale (D): ditto but use M-giant spectra if available
      if (off(2).eq.-1) then
         call bisect2(zlt,logts,nlogt,ix1,ix2)
         zlg = max(ming(ix1),ming(ix2))
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
         usegiant = testInMgiant(zlt)
c Use M-giant spectra alone if we are too far off in logg else interpolate to remain smooth
         if (zlg-logg.gt.0.5) then
            fmix = 1.d0
         else
            fmix = (zlg-logg)/0.5d0
         end if
         status = status + 8
         go to 100
      end if

c Off top of log t scale (E): use hottest model
      if (off(1).eq.1) then
         zlt = logts(nlogt)
         call bisect2(zlt,logts,nlogt,ix1,ix2)
         zlg = ming(ix1)
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
         status = status + 16
         go to 100
      end if

c Off bottom of log t scale (F): Interpolate in log g and hope for the best
      if (off(1).eq.-1) then
         zlg = logg
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
         zlt = max(mint(iy1),mint(iy2))
         call bisect2(zlt,logts,nlogt,ix1,ix2)
         status = status + 32
         go to 100
      end if

c At this point we are inside the bounds of the grid but we may not be out of the woods yet
c Now we know that we are comfortably within the grid, find the points to span
      zlt = logt
      zlg = logg
      call bisect2(zlt,logts,nlogt,ix1,ix2)
      call bisect2(zlg,loggs,nlogg,iy1,iy2)

c See whether we are safely inside the grid or not
      do i = 1, 2
         do j = 1, 2
            if (grid(1,ix1+i-1,iy1+j-1,iz1).gt.-1.d19) then
               tag(i,j) = 1
            else
               tag(i,j) = 0
            end if
         end do
      end do
      ttot = tag(1,1)+tag(2,1)+tag(1,2)+tag(2,2)
c If all four points are in the grid then great
      if (ttot.eq.4) then
         go to 100
      end if

c Check for being off the hot end of the grid (region (B))
      if (zlt.gt.maxt(iy1).or.zlt.gt.maxt(iy2)) then
c Find the best available log g
         zlg = max(ming(ix1),ming(ix2))
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
         status = status + 64
         go to 100
      else
c We are off the cold end of the grid :-(
         zlg = max(ming(ix1),ming(ix2))
c If we are close to cold dwarf models use those
         if (zlg-logg.lt.1.1d0) then
            call bisect2(zlg,loggs,nlogg,iy1,iy2)
            status = status + 128
            go to 100
         end if
c Cold giant models 
         zlg = logg
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
         zlt = max(mint(iy1),mint(iy2))
         call bisect2(zlt,logts,nlogt,ix1,ix2)
c If we are inside the M-giant grid then we're quite possibly better off using that
         usegiant = testInMgiant(logt)
         fmix = (1.d1**zlt - 1.d1**logt) / 2.5d2   ! (zlt-logt)/0.03d0   ! 0.03 === 250 K at that T
         fmix = max(0.d0,min(1.d0,fmix))
         status = status + 256
c         write (*,*) 'Usegiant=', usegiant, 'fmix=', fmix
         go to 100
      end if

c We should never get here
      write (*,*) 'CONFUSED in gt2mag: ', logg, logt, zlg,zlt,ix1,ix2,iy1,iy2
      STOP

100   continue
c We have now worked out where we are; calculate the f_i s
      fx = (zlt-logts(ix1))/(logts(ix2)-logts(ix1))
      fy = (zlg-loggs(iy1))/(loggs(iy2)-loggs(iy1))
      call interp3(fx,fy,fz,ix1,iy1,iz1,mags)
      if (.false.) then
         write (*,*) 'DEBUG', ix1, ix2, zlt, fx, iy1, iy2, zlg, fy, fz, (mags(i),i=1,nbands)
         if (mags(1).lt.-1.d10) then
            write (*,*) grid(1,ix1,iy1,iz1), grid(1,ix2,iy1,iz1)
            write (*,*) grid(1,ix1,iy2,iz1), grid(1,ix2,iy2,iz1)
            stop
         end if
      end if

c Sort out the crap with M-giants
      if (usegiant.ne.0) then
         call giantInterp(logt,zlm,magtmp,status)
         if (fmix.lt.0.d0 .or. fmix.gt.1.d0) then
            write (*,*) 'fmix ERROR', logg, logt, fmix
            stop
         end if
         do i = 1, nbands
            mags(i) = (1.d0-fmix) * mags(i) + fmix * magtmp(i)
         end do
            
c         if (usegiant.gt.1) then
c            do i = 1, nbands
c               mags(i) = magtmp(i)
c            end do
c         else
c            fz = logg/(-0.5d0)
c            fy = 1.d0 - fz
c            do i = 1, nbands
c               mags(i) = mags(i) * fy + magtmp(i) * fz
c            end do
c            status = status + 1024
c         end if
      end if
         

      return
      end

c Special arrangements to do giants correctly
      subroutine selectInterpolationGiant(logg, logt,
     :         zlg, zlt, ix1, ix2, iy1, iy2, debug)
      implicit none
      real*8 logg, logt, zlg, zlt
      integer ix1,ix2,iy1,iy2
      logical debug
      INCLUDE "ubvgrid.h"

      real*8 lowestPosblT,highestPosblT,epsT
      PARAMETER (epsT=1.d-5)    ! A small interval in logT

      debug = .false.
c Is our logg in the grid at all?
      zlg = logg
      if (logg.lt.loggs(1)) then
         iy1 = 1
         iy2 = 2
         debug = .true.
      else
c In the logg grid: find the correct place by interpolation
         call bisect2(zlg,loggs,nlogg,iy1,iy2)
      end if

c Now find the nearest temperature points for these loggs
      zlt = logt
      lowestPosblT = max(mint(iy1),mint(iy2)) + epsT
      highestPosblT = min(maxt(iy1),maxt(iy2)) - epsT
      if (zlt.lt.lowestPosblT) then
         call bisect2(lowestPosblT,logts,nlogt,ix1,ix2)
         debug = .true.
c         write (*,*) 'COLD', zlt, lowestPosblT, mint(iy1), mint(iy2), logts(ix1), logts(ix2)
      else if (zlt.gt.highestPosblT) then
         call bisect2(highestPosblT,logts,nlogt,ix1,ix2)
         debug = .true.
      else
         call bisect2(zlt,logts,nlogt,ix1,ix2)
      end if

c      write (*,*) 'GIANT', ix1, ix2, iy1, iy2

c Note that we do not change zlg,zlt so we may be extrapolating
      return
      end

c Test to see whether our temperature is inside the M-giant range
      INTEGER FUNCTION testInMgiant(zlt)
      implicit none
      include "ubvgrid.h"
      real*8 zlt

      testInMgiant = 0
      if (zlt.gt.giantTs(1).and.zlt.lt.giantTs(ngiant)) then
         testInMgiant = 1
      end if
      return
      end


c Search over the whole grid for the nearest model
      subroutine findNearestPoint(zlt, zlg, fz, iz1, iz2, mags)
      implicit none
      INCLUDE "ubvgrid.h"
      real*8 zlt,zlg,fz,mags(*)
      integer iz1, iz2

      real*8 d2min, d2
      integer i,j,l,ix,iy
      
      d2min = 1.e10
      do j = 1, nlogg
         do i = 1, nlogt
            if (grid(1,i,j,iz1).gt.-1.e19) then
               d2 = (zlt-logts(i))**2 + (zlg-loggs(j))**2
               if (d2.lt.d2min) then
                  d2min = d2
                  ix = i
                  iy = j
               end if
            end if
         end do
      end do
c      write (*,*) ix, iy, grid(1,ix,iy,iz1)
      do l = 1, nbands
         mags(l) = (1.d0-fz) * grid(l,ix,iy,iz1) + fz * grid(l,ix,iy,iz2)
      end do
      return
      end

c Interpolate in a cube of 8 points at i{x,y,z}1 + {0,1} 
      subroutine interp3(fx,fy,fz,ix1,iy1,iz1,mags)
      implicit none
      INCLUDE "ubvgrid.h"
      real*8 fx,fy,fz,mags(*)
      integer ix1,iy1,iz1

      integer i,j,k,l
      do i = 0, 1
         fz = 1.d0 - fz
         do j = 0, 1
            fy = 1.d0 - fy
            do k = 0, 1
               fx = 1.d0 - fx
               do l = 1, nbands
                  mags(l) = mags(l) + fx*fy*fz*
     :                      grid(l,ix1+k,iy1+j,iz1+i)
               end do
            end do
         end do
      end do
      return
      end


      subroutine hotInterp(logt,zlm,mags,status)
      implicit none
      INCLUDE "ubvgrid.h"
      real*8 logt,zlm,mags(nbands)
      integer status

      call interp1d(logt,zlm,mags,hotgrid,hotts,nhotmax,nhot,status,64)
      return
      end

      subroutine giantInterp(logt,zlm,mags,status)
      implicit none
      INCLUDE "ubvgrid.h"
      real*8 logt,zlm,mags(nbands)
      integer status

      integer i

c      write (*,*) 'giant interp: ngiant=', ngiant
c      write (*,*) (giantgrid(1,i,1),i=1,ngiant)

      call interp1d(logt,zlm,mags,giantgrid,giantts,ngiantmax,ngiant,status,512)
      return
      end


c Perform a general 1d interpolation (for either hot or giant grid)
      subroutine interp1d(x,zlm,mags,ygrid,yidx,nymax,ny,status,errnum)
      implicit none
      INCLUDE "ubvgrid.h"
      integer nymax,ny,status,errnum
      real*8 x,zlm,mags(nbands),ygrid(nbands,nymax,nmh),yidx(nymax)

      integer ix1,ix2,iz1,iz2,i,k,l
      real*8 fx,fz

c Check end
      if (x.lt.yidx(1)) then
         write (*,*) 'Off start of grid in interp1d', x,yidx(1)
         status = status + errnum
         x = yidx(1)
      else if (x.gt.yidx(ny)) then
         write (*,*) 'Off end of grid in interp1d', x,yidx(ny)
         status = status + errnum
         x = yidx(ny)
      end if
c      write (*,*) ny, nymax, (yidx(i),i=1,ny)
c      stop
      call bisect2(x,yidx,ny,ix1,ix2)

c Locate in metallicity
      call bisect2(zlm,mhs,nmh,iz1,iz2)

c Do the dance
      fx = (x-yidx(ix1))/(yidx(ix2)-yidx(ix1))
      fz = (zlm-mhs(iz1))/(mhs(iz2)-mhs(iz1))

      do l = 1, nbands
         mags(l) = 0.d0
      end do

      do i = 0, 1
         fz = 1.d0 - fz
         do k = 0, 1
            fx = 1.d0 - fx
            do l = 1, nbands
               mags(l) = mags(l) + fx*fz*
     :                   ygrid(l,ix1+k,iz1+i)
            end do
         end do
      end do
      return
      end

                     
         
c Find the indices i,j such that x(i)<y<x(j) where j-i=1 and x is a sorted
c array of nx doubles.
      subroutine bisect2(y,x,nx,i,j)
      implicit none
      real*8 y,x(*)
      integer nx,i,j

      integer k

      i=1
      j=nx
      do while (j-i.gt.1)
         k = (i+j)/2
         if (x(k).lt.y) then
            i=k
         else
            j=k
         end if
      end do
      return
      end
         





