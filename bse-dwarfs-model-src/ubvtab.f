***
      SUBROUTINE ubvtab
      implicit none
*
      INTEGER i,j,k,n
      INTEGER nzgr,ntgr,nggr
      PARAMETER(nzgr=8,ntgr=61,nggr= 11)
      INTEGER ntgr2,nggr2
      PARAMETER(ntgr2=91,nggr2=5)
*
      REAL feh
      REAL*8 zgr(nzgr),tgr(ntgr),ggr(nggr),ubv(nzgr,ntgr,nggr,5)
      COMMON /ubvdata/ zgr,tgr,ggr,ubv
      REAL*8 wtgr(ntgr2),wggr(nggr2),wubv(ntgr2,nggr2,5)
      COMMON /wubvdata/ wtgr,wggr,wubv
*
      OPEN(21,file='Kurucz.dat',
     &        form='formatted',status='old')
      do k = 1, nzgr
         do i = 1, ntgr
            do j = 1, nggr
               read(21,*)feh,tgr(i),ggr(j),(ubv(k,i,j,n),n=1,5)
            end do
            tgr(i) = log10(tgr(i))
         end do
c....... zgr=log(Z/0.02), assuming X=0.76-3*Z and Z(sun)=0.02
         zgr(k) = -log10((3.d0 + 37.425d0*10.d0**(-feh))/38.d0)
*        zgr(k) = -log10(0.07895 + 0.92105*10.0**(-feh))
      end do
      CLOSE(21)
*
      OPEN(22,file='wdhyd.dat',
     &        form='formatted',status='old')
      do j = 1,nggr2
         do i = 1,ntgr2
            read(22,*)wtgr(i),wggr(j),(wubv(i,j,k),k=1,5)
         enddo
      enddo
      do i = 1,ntgr2
         wtgr(i) = log10(wtgr(i))
      enddo
      CLOSE(22)
*
      RETURN
      END
***
