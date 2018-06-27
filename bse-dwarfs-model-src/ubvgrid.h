      INTEGER NlinesMax,NBands,Nmh,NlogTMax,NlogGMax,NHotMax,NGiantMax
c Lines in the input file
      PARAMETER (NlinesMax=10000)
c Dimensions of main grid
      PARAMETER (NMH=7,NlogTMax=200,NloggMax=20,NBands=3)
c Extra bits for hot stars and M giants (1d grids)
      PARAMETER (NHotMax=100,NgiantMax=100)

c Grid
      INTEGER nlogg,nlogt,jTminGiant
      REAL*8 loggs(NloggMax),logts(Nlogtmax),MHs(nMH)
      real*8 grid(nbands,nlogtmax,Nloggmax,nmh)
      real*8 mint(nlogtmax),maxt(nlogtmax)
      real*8 ming(nlogtmax),maxg(nlogtmax)
      real*8 logTminGiant

c Extra bits
      integer nhot,ngiant
      real*8 hotTs(nhotmax), hotGrid(nbands,nhotmax,nmh)
      real*8 giantTs(ngiantmax),giantGrid(nbands,ngiantmax,nmh)

      COMMON /ubvdata/ loggs,logts,mhs,grid,hotts,hotgrid,giantts,giantgrid,logTminGiant,mint,maxt,ming,maxg,
     :                 nlogg,nlogt,nhot,ngiant,jTminGiant
