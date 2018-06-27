      PROGRAM testbcs
      IMPLICIT NONE
      INCLUDE "ubvgrid.h"
      integer i,k,status
      real*8 t,logr,logt,logl,m,logg,mags(nbands),magsun

      integer ix1,iy1,iz1,j
      common /spy/ ix1,iy1,iz1
      
      call loadtab

      open (10,file="plot.1")
      open (11,file="mags.out")
      do
         read (10,*,end=100) k, t, logr, logt, logl, m
         logg = log10(m) - 2.d0*logr + 4.4383174
         magsun = 4.75d0 - 2.5d0*logl
         call gt2mag(logg,logt,0.d0,mags,status)
         write (11,99001) k, t, logr, logt, logl, m, logg,magsun,(magsun-mags(i),i=1,nbands),status,ix1,iy1,iz1
         call flush(11)
      end do
100   close(10)
      close(11)

      open (11,file="grid.out")
      do i = 1, nlogg
         write (11,*) (grid(1,j,i,1),j=1,nlogt)
      end do
      close(11)

99001 format(i5,1x,1e15.6,14(1x,f10.6),4(1x,i3))
      end


