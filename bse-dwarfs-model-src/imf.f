***
      real*8 FUNCTION IMF1(X)
      IMPLICIT NONE
      real*8 X,G1,G2,G3,G4
      DATA  G1,G2,G3,G4 /0.19,1.55,0.050,0.6/
*
* A function to generate masses from the IMF of Kroupa, Tout
* & Gilmore, MNRAS, 1993, 262, 545.
*
      IMF1 = 0.08 + (G1*X**G2 + G3*X**G4)/(1.0 - X)**0.58
*
      RETURN
      END
***
      real*8 FUNCTION IMF2(X,M0,ALPHA)
      IMPLICIT NONE
      real*8 X,M0,ALPHA
*
* A function to generate masses from a Salpeter IMF.
*
      IMF2 = M0/(1.0 - X)**ALPHA
*
      RETURN
      END
***
      real*8 FUNCTION IMF3(X)
      IMPLICIT NONE
      real*8 X,XX,G1,G2,G3,G4
      DATA  G1,G2,G3,G4 /0.75,0.04,0.25,1.04/
*
* A function to generate masses from the IMF of Kroupa, Tout
* & Gilmore, MNRAS, 1991, 251, 293.
*
      XX = 1.0 - X
      IMF3 = 0.33*((1.0/(XX**G1 + G2*XX**G3)) - (XX**2/G4))
*
      RETURN
      END
***
