c Power law mass ratio
 
      SUBROUTINE getm2(m1,m2,tb)
      IMPLICIT NONE
      REAL*8 m1,m2,gammaq,gammaql,gammaqs,qb,X,Xb,a1,tb

      include "rand.h"
      REAL*8 ran3

      qb = 0.3d0
      X = ran3(seed)

      if (m1.lt.0.1d0) then
            gammaq = 4.2d0
            a1 = (1.d0 + gammaq)/(1.d0 - 0.1d0**(1.d0 + gammaq))
c            m2 = m1*X**(1.d0/(1.d0+gammaq))
            m2 = m1*(-(((1.d0 + gammaq)*(-((0.1d0**(1.d0 + gammaq)*a1)
     &      /(1.d0 + gammaq)) - X))/a1))**(1.d0/(1.d0 + gammaq))
      else if (m1.ge.0.1d0 .and. m1.lt.0.3d0) then
            gammaq = 1.9d0
            a1 = (1.d0 + gammaq)/(1.d0 - 0.1d0**(1.d0 + gammaq))
c            m2 = m1*X**(1.d0/(1.d0+gammaq))
            m2 = m1*(-(((1.d0 + gammaq)*(-((0.1d0**(1.d0 + gammaq)*a1)
     &      /(1.d0 + gammaq)) - X))/a1))**(1.d0/(1.d0 + gammaq))
      else if (m1.ge.0.3d0 .and. m1.lt.0.8) then
            gammaq = -0.2d0
            a1 = (1.d0 + gammaq)/(1.d0 - 0.1d0**(1.d0 + gammaq))
c            m2 = m1*X**(1.d0/(1.d0+gammaq))
            m2 = m1*(-(((1.d0 + gammaq)*(-((0.1d0**(1.d0 + gammaq)*a1)
     &      /(1.d0 + gammaq)) - X))/a1))**(1.d0/(1.d0 + gammaq))



      else if (m1.ge.0.8d0 .and. m1.lt.1.6d0) then
            gammaqs = 0.3d0
            if (log10(tb).le.6.d0) then
                  gammaql = -0.5d0
                  Xb = 0.26151d0
            else if (log10(tb).gt.6.d0) then
                  gammaql = -1.1d0
                  Xb = 0.34018d0
            end if

           a1 = -1.d0/((-qb**(-gammaql + gammaqs)*(1.d0 - qb**(1.d0 
     &       + gammaql)))/(1.d0 + gammaql) - ((-0.1d0**(1.d0 + 
     &        gammaqs) + qb**(1.d0 + gammaqs)))/(1.d0 + gammaqs))

            if (X.le.Xb) then
                  m2 = m1*(((1.d0 + gammaqs)*(X + (0.1d0**(1.d0 + 
     &             gammaqs)*a1)/(1.d0 + gammaqs)))/a1)**(1.d0/(1.d0 + 
     &              gammaqs))
            else
                  m2 = m1*((qb**(gammaql - gammaqs)*(1.d0 + 
     &             gammaql)*(-1.d0 + X + (a1*qb**(-gammaql + 
     &              gammaqs))/(1.d0 + gammaql)))/a1)**(1.d0/(1.d0 
     &               + gammaql))
            end if

      else if (m1.ge.1.6d0 .and. m1.lt.5.d0) then
            if (log10(tb).le.2.d0) then
                  gammaql = -0.5d0
                  gammaqs =  0.2d0
                  Xb = 0.269847d0
            else if (log10(tb).gt.2.d0 .and. log10(tb).le.4.d0) then
                  gammaql = -0.9d0
                  gammaqs =  0.1d0
                  Xb = 0.332591d0
            else if (log10(tb).gt.4.d0 .and. log10(tb).le.6.d0) then
                  gammaql = -1.4d0
                  gammaqs = -0.5d0
                  Xb = 0.106698d0
            else if (log10(tb).gt.6.d0) then    
                  gammaql = -2.0d0
                  gammaqs = -1.0d0
                  Xb = 0.278053d0
            end if
!            1.d0/((qb**(-1.d0 - gammaql)*(1.d0 - qb**(1.d0 + gammaql)))/(1.d0 + gammaql) - Log(0.1d0) + Log(qb))


            if (gammaqs.eq.-1.0d0) then
                  a1 = 1.d0/((qb**(-1.d0 - gammaql)*(1.d0 - qb**(1.d0
     &              + gammaql)))/(1.d0 + gammaql)
     &               - Log(0.1d0) + Log(qb))
            else
                  a1 = -1.d0/((-qb**(-gammaql + gammaqs)*(1.d0 
     &              - qb**(1.d0 + gammaql)))/(1.d0 + gammaql) - 
     &             ((-0.1d0**(1.d0 + gammaqs) +
     &             qb**(1.d0 + gammaqs)))/(1.d0 + gammaqs))
            end if

            if (X.le.Xb) then
                  if (gammaqs.eq.-1.0d0) then
                        m2 = m1*0.1d0*Exp(X/a1)
                  else
                        m2 = m1*(((1.d0 + gammaqs)*(X + (0.1d0**(1.d0 
     &                   + gammaqs)*a1)/(1.d0 + 
     &                    gammaqs)))/a1)**(1.d0/(1.d0 + gammaqs))
                  end if
            else
                  m2 = m1*((qb**(gammaql - gammaqs)*(1.d0 + 
     &             gammaql)*(-1.d0 + X + (a1*qb**(-gammaql + 
     &              gammaqs))/(1.d0 + gammaql)))/a1)**(1.d0/(1.d0 
     &               + gammaql))
            end if

      else if (m1.ge.5.d0) then
            if (log10(tb).le.2.d0) then
                  gammaql = -0.5d0
                  gammaqs =  0.1d0
                  Xb = 0.278536d0
            else if (log10(tb).gt.2.d0 .and. log10(tb).le.4.d0) then
                  gammaql = -1.7d0
                  gammaqs = -0.2d0
                  Xb = 0.473257d0
            else if (log10(tb).gt.4.d0 .and. log10(tb).le.6.d0) then
                  gammaql = -2.0d0
                  gammaqs = -1.2d0
                  Xb = 0.637053d0
            else if (log10(tb).gt.6.d0) then
                  gammaql = -2.0d0
                  gammaqs = -1.5d0
                  Xb = 0.67654d0
            end if

            a1 = -1.d0/((-qb**(-gammaql + gammaqs)*(1.d0 - qb**(1.d0 
     &       + gammaql)))/(1.d0 + gammaql) - ((-0.1d0**(1.d0 + 
     &        gammaqs) + qb**(1.d0 + gammaqs)))/(1.d0 + gammaqs))

            if (X.le.Xb) then
                  m2 = m1*(((1.d0 + gammaqs)*(X + (0.1d0**(1.d0 + 
     &             gammaqs)*a1)/(1.d0 + gammaqs)))/a1)**(1.d0/(1.d0  
     &             + gammaqs))
            else
                  m2 = m1*((qb**(gammaql - gammaqs)*(1.d0 + 
     &             gammaql)*(-1.d0 + X + (a1*qb**(-gammaql + 
     &              gammaqs))/(1.d0 + gammaql)))/a1)**(1.d0/(1.d0 
     &               + gammaql))
            end if
      end if


      RETURN
      END
