      FUNCTION lb2red(l, b)
      IMPLICIT NONE

      REAL*8 l, b, x1, x2, y1, y2, lb2red, readRed, interpRed
      REAL*8 ROUND
      INTEGER idl_base, idb_base, idx, i, getIndex

      PARAMETER (ROUND=0.5d0)

c Calculate x values for interpolation
      x1 = INT((l + (ROUND/2.0d0)) * (1d0/ROUND)) / (1d0/ROUND)
      if (x1.le.l .AND. x1.gt.0d0) then
        x2 = x1 - ROUND
      else if (x1.ge.l .AND. x1.lt.360d0) then
        x2 = x1 + ROUND
      else if (x1.eq.0d0) then
        x2 = ROUND
      else if (x1.eq.360d0) then
        x2 = 360d0 - ROUND
      end if

c Calculate y values for interpoloation
      y1 = INT((b + (ROUND/2.0d0)) * (1d0/ROUND)) / (1d0/ROUND)
      if (y1.le.b .AND. y1.gt.0d0) then
        y2 = y1 - ROUND
      else if (y1.ge.b .AND. y1.lt.180d0) then
        y2 = y1 + ROUND
      else if (y1.eq.0d0) then
        y2 = ROUND
      else if (x1.eq.180d0) then
        y2 = 180d0 - ROUND
      end if

      lb2red = interpRed(x1, x2, y1, y2, l, b)
      
      END FUNCTION lb2red

c Function to read reddening values from file
      FUNCTION readRed(x, y)
      
      REAL*8 bad_value, readRed, x, y
      INTEGER idx, getIndex

      idx = getIndex(x, y)

      OPEN(unit=120, file='reddening.in', action='read')
      
      do i=1, idx
        read(120,*) readRed
      end do
      close(120)

      END FUNCTION readred

c Function to get index of rounded l and b values. Assumes resoluation of 0.5deg
      FUNCTION getIndex(x, y)

      REAL*8 x, y, ROUND
      INTEGER idl_base, idb_base, getIndex

      PARAMETER (ROUND=0.5d0)

      idl_base = (x * (1d0/ROUND)) * (180d0/ROUND + 1) 
      idb_base = (y * (1d0/ROUND)) + 1
      getIndex = idl_base + idb_base

      END FUNCTION getIndex

c Function to interpolate between the 4 points
      FUNCTION interpRed(x1, x2, y1, y2, l, b)

      REAL*8 interpRed, x1, x2, y1, y2, l, b
      REAL*8 val11, val12, val21, val22, vale
      REAL*8 readRed
      INTEGER getIndex

      val11 = readRed(x1, y1)
      val12 = readRed(x1, y2)
      val21 = readRed(x2, y1)
      val22 = readRed(x2, y2)

      interpRed = (1d0 / ((x2 - x1)*(y2 - y1))) * 
     +       ((val11 * (x2 - l) * (y2 - b)) + 
     +       (val21 * (l - x1) * (y2 - b)) + 
     +       (val12 * (x2 - l) * (b - y1)) + 
     +       (val22 * (l - x1) * (b - y1)))


      END FUNCTION interpRed