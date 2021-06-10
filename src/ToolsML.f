c=======================================================================
c=======================================================================
c     TOOLS: (M)ARGINAL (L)IKELIHOODS
c
c     SUBROUTINE NAMES:
c     · LOGNORNIG
c     · LOGPOIGAM
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine lognornig(nobs,obs,npars,pars,labels,indj,val)
c=======================================================================
c=======================================================================
c     BEGIN: LOGNORNIG SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TIME SERIE (Y[T]:T=1,...,N), A SET OF CLUSTER LABELS
c     (E[T]:T=1,...,N) WITH 1=E[1]<=E[2]<=···<=E[N], AND
c     J IN {1,...,E[N]}, "LOGNORNIG" RETURNS THE LOG MARGINAL
c     LIKELIHOOD ML(Y[T]:T IN S[J]|THETA) INDUCED BY THE FOLLOWING
c     MODEL:
c
c     A) Y[T]|M,S^2~NORMAL(M,S^2):T IN S[J]={I IN {1,...,N}:E[I]=J}.
c
c        · GIVEN M AND S^2, OBSERVATIONS ARE INDEPENDENT.
c
c     B) M|S^2~NORMAL(MU,{(S^2)/KAPPA}).
c
c        · MEAN MU AND PRECISION KAPPA ARE FIXED.
c
c     C) S^2~INVERSE-GAMMA(ALPHA,BETA).
c
c        · SHAPE ALPHA AND SCALE BETA ARE FIXED.
c
c     IN THIS CASE, THETA=(MU,KAPPA,ALPHA,BETA). LETING D=CARD(S[J]),
c     ML(Y[T]:T IN S[J]|THETA)} IS THE DENSITY OF A D-DIMENSIONAL
c     T DISTRIBUTION WITH DEGREES OF FREEDOM 2*ALPHA, LOCATION VECTOR
c     MU*1[D] AND SCALE MATRIX (BETA/ALPHA)*(I[D]+(1/KAPPA)*1[D]1[D]').
c
c     1) 1[D]: D-DIMENSIONAL VECTOR WITH ALL ENTRIES EQUAL 1.
c
c     2) I[D]: D-DIMENSIONAL IDENTITY MATRIX.
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TIME SERIES (Y[T]:T=1,...,N)
      integer nobs
c     obs: TIME SERIES (Y[T]:T=1,...,N)
      real*8 obs(nobs)
c     npars: MAXIMUM NUMBER OF ADMISSIBLE PARAMETERS
      integer npars
c     pars: PARAMETER THETA=(MU,KAPPA,ALPHA,BETA)
      real*8 pars(npars)
c     labels: CLUSTER LABELS (E[T]:T=1,...,N)
      integer labels(nobs)
c     indj: INDEX J
      integer indj
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: LOG{ML(Y[T]:T IN S[J]|THETA)}
      real*8 val
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer s
c     EXCLUSIVE FOR CONSTANT LOG(PI)
      real*8 logpi
c     EXCLUSIVE FOR STORING D
      integer d
c     EXCLUSIVE FOR STORING LOG{ML(Y[T]:T IN S[J]|THETA)}
      real*8 mu
      real*8 kappa
      real*8 alpha
      real*8 beta
      real*8 sy
      real*8 s2y
c     OTHERS
      real*8 dw
      real*8 rw
      real*8 sw
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      if (npars.lt.4) then
         print *,'Error in lognornig subroutine: check npars'
         stop
      end if
      logpi=1.144729885849400174143427351353d0
      d=0
      sy=0.d0
      s2y=0.d0
      do s=1,nobs
         if (labels(s).eq.indj) then
            d=d+1
            sy=sy+obs(s)
            s2y=s2y+(obs(s)*obs(s))
         end if
      end do
      mu=pars(1)
      kappa=pars(2)
      alpha=pars(3)
      beta=pars(4)
c=======================================================================
c     ALGORITHM
c=======================================================================
      dw=dble(d)
      val=dlgama(alpha+(0.5d0*dw))-dlgama(alpha)
      rw=0.5d0*(dlog(kappa)-dlog(kappa+dw))
      sw=(dw*0.5d0)*(logpi+dlog(2.d0*beta))
      val=val+(rw-sw)
      rw=((mu*mu)*kappa)+s2y
      sw=(((kappa*mu)+sy)*((kappa*mu)+sy))/(kappa+dw)
      val=val-(alpha+(0.5d0*dw))*dlog(1.d0+((0.5d0*(rw-sw))/beta))
c=======================================================================
      return
c     END: LOGNORNIG SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logpoigam(nobs,obs,npars,pars,labels,indj,val)
c=======================================================================
c=======================================================================
c     BEGIN: LOGPOIGAM SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TIME SERIE (Y[T]:T=1,...,N), A SET OF CLUSTER LABELS
c     (E[T]:T=1,...,N) WITH 1=E[1]<=E[2]<=···<=E[N], AND
c     J IN {1,...,E[N]}, "LOGNORNIG" RETURNS THE LOG MARGINAL
c     LIKELIHOOD ML(Y[T]:T IN S[J]|THETA) INDUCED BY THE FOLLOWING
c     MODEL:
c
c     A) Y[T]|L~POISSON(L):T S[J]={I IN {1,...,N}:E[I]=J}.
c
c        · GIVEN L, OBSERVATIONS ARE INDEPENDENT.
c
c     B) L~GAMMA(ALPHA,BETA).
c
c        · SHAPE ALPHA AND RATE BETA ARE FIXED.
c
c     IN THIS CASE, THETA=(ALPHA,BETA). LETING D=CARD(S[J]),
c     ML(Y[T]:T IN S[J]|THETA)} IS THE DENSITY OF A NEGATIVE
c     MULTINOMIAL DISTRIBUTION WITH DISPERSION ALPHA AND
c     PROBABILITY EVENTS {1/(B+N)}*1[D].
c
c     1) 1[D]: D-DIMENSIONAL VECTOR WITH ALL ENTRIES EQUAL 1.
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TIME SERIES (Y[T]:T=1,...,N)
      integer nobs
c     obs: TIME SERIES (Y[T]:T=1,...,N)
      real*8 obs(nobs)
c     npars: MAXIMUM NUMBER OF ADMISSIBLE PARAMETERS
      integer npars
c     pars: PARAMETER THETA=(MU,KAPPA,ALPHA,BETA)
      real*8 pars(npars)
c     labels: CLUSTER LABELS (E[T]:T=1,...,N)
      integer labels(nobs)
c     indj: INDEX J
      integer indj
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: LOG{ML(Y[T]:T IN S[J]|THETA)}
      real*8 val
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer s
c     EXCLUSIVE FOR STORING D
      integer d
c     EXCLUSIVE FOR STORING LOG{ML(Y[T]:T IN S[J]|THETA)}
      real*8 alpha
      real*8 beta
      real*8 sy
      real*8 slgy
c     OTHERS
      real*8 dw
      real*8 rw
      real*8 sw
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      if (npars.lt.2) then
         print *,'Error in logpoigam subroutine: check npars'
         stop
      end if
      d=0
      sy=0.d0
      slgy=0.d0
      do s=1,nobs
         if (labels(s).eq.indj) then
            d=d+1
            sy=sy+obs(s)
            slgy=slgy+dlgama(obs(s)+1.d0)
         end if
      end do
      alpha=pars(1)
      beta=pars(2)
c=======================================================================
c     ALGORITHM
c=======================================================================
      dw=dble(d)
      rw=dlgama(alpha+sy)-(dlgama(alpha)+slgy)
      sw=(alpha*(dlog(beta)-dlog(beta+dw)))-(sy*dlog(beta+dw))
      val=rw+sw
c=======================================================================
      return
c     END: LOGPOIGAM SUBROUTINE
      end
c=======================================================================
c=======================================================================

