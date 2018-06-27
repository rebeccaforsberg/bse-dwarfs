***
      subroutine lt2ubv(logl,logt,mass,logz,bolc,mv,uminb,bminv,vmini)
      implicit none
c.... computes values of Mv, U-B, B-V and V-I for given log L, log T,
c.... mass and log(Z/0.02)
*
      integer k,ig,it,iz
      integer nzgr,ntgr,nggr
      parameter(nzgr=8,ntgr=61,nggr=11)
      integer indx
      external indx
      real*8 zgr(nzgr),tgr(ntgr),ggr(nggr),ubv(nzgr,ntgr,nggr,5)
      common /ubvdata/ zgr,tgr,ggr,ubv
      real*8 cm(5),gconst
      parameter(gconst=-10.6071d0)
      real*8 logl,logt,mass,logz,mv,uminb,bminv,vmini
      real*8 logm,logg,dg1,dg2,dt1,dt2,dz1,dz2,mbol,bolc
*
      logm = dlog10(mass)
      logg = logm + 4.d0*logt - logl + gconst
c.... find indices of log Z, log g and log T to interpolate between.
c.... don't allow extrapolation outside log Z and log g grid.
      ig = indx(logg,ggr,nggr)
      it = indx(logt,tgr,ntgr)
      iz = indx(logz,zgr,nzgr)
      dg1 = (logg - ggr(ig-1))/(ggr(ig) - ggr(ig-1))
      dg1 = max(0.d0, min(1.d0, dg1))
      dg2 = 1.d0 - dg1
      dt1 = (logt - tgr(it-1))/(tgr(it) - tgr(it-1))
      dt2 = 1.d0 - dt1
      dz1 = (logz - zgr(iz-1))/(zgr(iz) - zgr(iz-1))
      dz1 = max(0.d0, min(1.d0, dz1))
      dz2 = 1.d0 - dz1
      do k = 1, 5
         cm(k) = ((ubv(iz,it,ig,k)*dg1 + ubv(iz,it,ig-1,k)*dg2)*dt1
     &            + (ubv(iz,it-1,ig,k)*dg1 +
     &            ubv(iz,it-1,ig-1,k)*dg2)*dt2)*dz1 +
     &           ((ubv(iz-1,it,ig,k)*dg1 +
     &            ubv(iz-1,it,ig-1,k)*dg2)*dt1 +
     &           (ubv(iz-1,it-1,ig,k)*dg1 +
     &            ubv(iz-1,it-1,ig-1,k)*dg2)*dt2)*dz2
      enddo
      mbol = 4.75d0 - 2.5d0*logl
      bolc = cm(1)
      mv = mbol - bolc
      uminb = cm(2)
      bminv = cm(3)
      vmini = cm(4) + cm(5)
*
      return
      end
***
      integer function indx(ax,xx,nx)
c.....finds index of ax in monotonously increasing or decreasing array xx
      implicit none
      integer nx,j,jl,jh
      real*8 ax,xx(nx),sx
*
      sx = xx(nx) - xx(1)
      jl = 1
      jh = nx
 1    if (jh-jl.gt.1) then
         j = (jh + jl)/2
         if ((ax-xx(j))*sx.gt.0.d0) then
            jl = j
         else
            jh = j
         end if
         goto 1
      end if
      indx = jh
*
      return
      end
***
