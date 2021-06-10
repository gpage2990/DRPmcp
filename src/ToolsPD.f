c=======================================================================
c=======================================================================
c     TOOLS: (P)ROBABILITY (D)ISTRIBUTIONS
c
c     SUBROUTINE NAMES:
c     · MVTLOGITD
c     · MVTLOGITR
c     · LOGMVTD
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logmvtd(dmn,x,nu,mu,invsigma,val)
c=======================================================================
c=======================================================================
c     BEGIN: LOGMVTD SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "LOGMVTD" RETURNS THE UNNORMALIZED LOG DENSITY, EVALUATED AT
c     X IN (0,1)^D, OF A D-DIMENSIONAL T DISTRIBUTION WITH DEGREES OF
c     FREEDOM NU, LOCATION VECTOR MU AND SCALE MATRIX SIGMA.
c=======================================================================
c     INPUTS
c=======================================================================
c     dmn: DIMENSION OF THE MULTIVARIATE T
      integer dmn
c     x: VECTOR WHERE THE DENSITY WILL BE EVALUATED
      real*8 x(dmn)
c     nu: DEGREES OF FREEDOM NU
      real*8 nu
c     mu: LOCATION VECTOR MU
      real*8 mu(dmn)
c     invsigma: INVERSE OF SCALE MATRIX SIGMA
      real*8 invsigma(dmn,dmn)
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: LOG{T[D](X|NU,MU,SIGMA)}
      real*8 val
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
      integer jj
c     OTHERS
      real*8 rw
      real*8 sw
c=======================================================================
c     ALGORITHM
c=======================================================================
      sw=0.d0
      do ii=1,dmn
         do jj=1,dmn
            rw=((x(ii)-mu(ii))*invsigma(ii,jj))*(x(jj)-mu(jj))
            sw=sw+rw
         end do
      end do
      val=(0.5d0*(-nu-dble(dmn)))*dlog(1.d0+(sw/nu))
c=======================================================================
      return
c     END: LOGMVTD SUBROUTINE
      end
c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine mvtlogitd(dmn,x,nu,mu,invsigma,val)
c=======================================================================
c=======================================================================
c     BEGIN: MVTLOGITD SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "MVTLOGITD" RETURNS THE LOG DENSITY, EVALUATED AT X IN (0,1)^D
c     AND IGNORING THE NORMALIZING CONSTANT, OF A D-DIMENSIONAL LOGIT-T
c     DISTRIBUTION WITH DEGREES OF FREEDOM NU, LOCATION VECTOR MU AND
c     SCALE MATRIX SIGMA.
c=======================================================================
c     INPUTS
c=======================================================================
c     dmn: DIMENSION OF THE MULTIVARIATE LOGIT-T
      integer dmn
c     x: VECTOR WHERE THE DENSITY WILL BE EVALUATED
      real*8 x(dmn)
c     nu: DEGREES OF FREEDOM NU
      real*8 nu
c     mu: LOCATION VECTOR MU
      real*8 mu(dmn)
c     invsigma: INVERSE OF SCALE MATRIX SIGMA
      real*8 invsigma(dmn,dmn)
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: LOG{LOGIT-T[D](X|NU,MU,SIGMA)}
      real*8 val
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
      integer jj
c     OTHERS
      real*8 jw
      real*8 rw
      real*8 sw
      real*8 yw(dmn)
c=======================================================================
c     ALGORITHM
c=======================================================================
      jw=0.d0
      do ii=1,dmn
c        YW=(LOG{X[I]/(1-X[I])}:I=1,...,D)
         yw(ii)=dlog(x(ii))-dlog(1.d0-x(ii))
c        JW=SUM(LOG(X[I])+LOG(1-X[I]):I=1,...,D)
         jw=jw+(dlog(x(ii))+dlog(1.d0-x(ii)))
      end do
      sw=0.d0
      do ii=1,dmn
         do jj=1,dmn
            rw=((yw(ii)-mu(ii))*invsigma(ii,jj))*(yw(jj)-mu(jj))
            sw=sw+rw
         end do
      end do
      val=((0.5d0*(-nu-dble(dmn)))*dlog(1.d0+(sw/nu)))-jw
c=======================================================================
      return
c     END: MVTLOGITD SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine mvtlogitr(dmn,nu,mu,cholsigma,val)
c=======================================================================
c=======================================================================
c     BEGIN: MVTLOGITR SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "MVTLOGITR" RETURNS A RANDOM VECTOR GENERATED FROM A
c     D-DIMENSIONAL LOGIT-T DISTRIBUTION WITH DEGREES OF FREEDOM NU,
c     LOCATION VECTOR MU AND SCALE MATRIX SIGMA.
c=======================================================================
c     INPUTS
c=======================================================================
c     dmn: DIMENSION OF THE MULTIVARIATE LOGIT-T
      integer dmn
c     nu: DEGREES OF FREEDOM NU
      real*8 nu
c     mu: LOCATION VECTOR MU
      real*8 mu(dmn)
c     cholsigma: CHOLESKY DECOMPOSITION OF SCALE MATRIX SIGMA
      real*8 cholsigma(dmn,dmn)
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: X~LOGIT-T[D](NU,MU,SIGMA)
      real*8 val(dmn)
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
      real*8 chisqr
      real*8 normr
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
      integer jj
c     OTHERS
      real*8 rw
      real*8 sw
      real*8 tw(dmn)
      real*8 zw(dmn)
c=======================================================================
c     ALGORITHM
c=======================================================================
      sw=chisqr(nu)
      sw=dsqrt(nu/sw)
      do ii=1,dmn
         zw(ii)=normr(0.d0,1.d0)
      end do
      do ii=1,dmn
         rw=0.d0
         do jj=1,dmn
            rw=rw+(cholsigma(ii,jj)*zw(jj))
         end do
c        TW~T[D](NU,MU,SIGMA)
         tw(ii)=mu(ii)+(sw*rw)
      end do
      do ii=1,dmn
         zw(ii)=tw(ii)-dlog(1.d0+dexp(tw(ii)))
         val(ii)=dexp(zw(ii))
      end do
c=======================================================================
      return
c     END: MVTLOGITR SUBROUTINE
      end
c=======================================================================
c=======================================================================

