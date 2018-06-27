! generates number from gaussian distribution of mean 0 and variance 1
      FUNCTION gasdev()
      INTEGER idum
      REAL*8 gasdev
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran3
      include "rand.h"
      SAVE iset,gset
      DATA iset/0/
      IF (iset.eq.0) THEN
1       v1=2.*ran3(seed)-1.
        v2=2.*ran3(seed)-1.
        rsq=v1**2+v2**2
        IF(rsq.ge.1. .OR. rsq.eq.0.)GOTO 1
        fac=SQRT(-2.*LOG(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      ELSE
        gasdev=gset
        iset=0
      END IF
      RETURN
      END 





