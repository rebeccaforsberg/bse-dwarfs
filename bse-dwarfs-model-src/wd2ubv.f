***
      subroutine wd2ubv(logl,logt,mass,bolc,mv,uminb,bminv,vmini)
c.... computes values of Mv, U-B, B-V and V-I for given log L, log T,
c.... mass and log(Z/0.02)
      implicit none
      integer k,ig,it
      integer ntgr,nggr
      parameter(ntgr=91,nggr=5)
      integer indx
      external indx
      real*8 tgr(ntgr),ggr(nggr),ubv(ntgr,nggr,5)
      common /wubvdata/ tgr,ggr,ubv
      real*8 cm(5),gconst
      parameter(gconst = -10.6071d0)
      real*8 mbol,bolc,mv,uminb,bminv,vmini
      real*8 logm,logg,dg1,dg2,dt1,dt2,mass,logt,logl
*
      logm = dlog10(mass)
      logg = logm + 4.d0*logt - logl + gconst
c.... find indices of log g and log T to interpolate between.
c.... don't allow extrapolation outside log g grid.
      ig = indx(logg,ggr,nggr)
      it = indx(logt,tgr,ntgr)
      dg1 = (logg - ggr(ig-1))/(ggr(ig) - ggr(ig-1))
      dg1 = MAX(0.d0,MIN(1.d0,dg1))
      dg2 = 1.d0 - dg1
      dt1 = (logt - tgr(it-1))/(tgr(it) - tgr(it-1))
      dt2 = 1.d0 - dt1
      do k = 1,5
         cm(k) = (ubv(it,ig,k)*dg1 + ubv(it,ig-1,k)*dg2)*dt1
     &         + (ubv(it-1,ig,k)*dg1 + ubv(it-1,ig-1,k)*dg2)*dt2
      enddo
      mbol = 4.75d0 - 2.5d0*logl
      bolc = cm(1)
      mv = mbol - bolc
      uminb = cm(2)
      bminv = cm(3)
      vmini = cm(4) + cm(5)
      end
***
***
