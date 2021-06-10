c=======================================================================
c=======================================================================
c     TOOLS: (G)IBBS (S)AMPLING
c
c     SUBROUTINE NAMES:
c     · LOWERINDEX
c     · UPPERINDEX
c     · LOGML
c     · RWMHP
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine lowerindex(ncs,cs,indt,indat)
c=======================================================================
c=======================================================================
c     BEGIN: LOWERINDEX SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A BINARY SEQUENCE C=(C[S]:S=1,...,N-1) IN {0,1}^(N-1)
c     AND 0<T<N, "LOWERINDEX" RETURNS THE INDEX A[T] DEFINED AS
c
c     A[T]=
c
c     · MAX{S=1,...,T-1:C[S]=1}, IF T>1 AND
c       SUM(C[S]:S=1,...,T-1)>0.
c
c     · 0, OTHERWISE.
c=======================================================================
c     INPUTS
c=======================================================================
c     ncs: ARGUMENT N
      integer ncs
c     cs: ARGUMENT C=(C[S]:S=1,...,N-1)
      integer cs(ncs-1)
c     indt: ARGUMENT T
      integer indt
c=======================================================================
c     OUTPUTS
c=======================================================================
c     indat: INDEX A[T]
      integer indat
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ss
c     OTHERS
      integer sw
c=======================================================================
c     ALGORITHM
c=======================================================================
c     · CONDITION T=1 IMPLIES A[T]=0
      if (indt.eq.1) then
         indat=0
      else
c        SW=SUM(C[S]:S=1,...,T-1)
         sw=0
         do ss=1,(indt-1)
            sw=sw+cs(ss)
         end do
c        · CONDITION SUM(C[S]:S=1,...,T-1)=0 IMPLIES A[T]=0
         if (sw.eq.0) then
            indat=0
         else
c           A[T]=MAX{S=1,...,T-1:C[S]=1}
            indat=indt-1
            ss=indt-1
            do while (ss.gt.0)
               if (cs(ss).eq.0) then
                  indat=indat-1
                  ss=ss-1
               else
                  ss=0
               end if
            end do
         end if
      end if
c=======================================================================
      return
c     END: LOWERINDEX SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine upperindex(ncs,cs,indt,indbt)
c=======================================================================
c=======================================================================
c     BEGIN: UPPERINDEX SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A BINARY SEQUENCE C=(C[S]:S=1,...,N-1) IN {0,1}^(N-1)
c     AND 0<T<N, "UPPERINDEX" RETURNS THE INDEX B[T] DEFINED AS
c
c     B[T]=
c
c     · MIN{S=T+1,...,N-1:C[S]=1}, IF T<N-1 AND
c       SUM(C[S]:S=T+1,...,N-1)>0.
c
c     · N, OTHERWISE.
c=======================================================================
c     INPUTS
c=======================================================================
c     ncs: ARGUMENT N
      integer ncs
c     cs: ARGUMENT C=(C[S]:S=1,...,N-1)
      integer cs(ncs-1)
c     indt: ARGUMENT T
      integer indt
c=======================================================================
c     OUTPUTS
c=======================================================================
c     indbt: INDEX B[T]
      integer indbt
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ss
c     OTHERS
      integer sw
c=======================================================================
c     ALGORITHM
c=======================================================================
c     · CONDITION T=N-1 IMPLIES B[T]=N
      if (indt.eq.(ncs-1)) then
         indbt=ncs
      else
c        SW=SUM(C[S]:S=T+1,...,N-1)
         sw=0
         do ss=(indt+1),(ncs-1)
            sw=sw+cs(ss)
         end do
c        · CONDITION SUM(C[S]:S=T+1,...,N-1)=0 IMPLIES B[T]=N
         if (sw.eq.0) then
            indbt=ncs
         else
c           B[T]=MIN{S=T+1,...,N-1:C[S]=1}
            indbt=indt+1
            ss=indt+1
            do while (ss.gt.0)
               if (cs(ss).eq.0) then
                  indbt=indbt+1
                  ss=ss+1
               else
                  ss=0
               end if
            end do
         end if
      end if
c=======================================================================
      return
c     END: UPPERINDEX SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logml(indtype,nobs,obs,npars,pars,labels,indj,val)
c=======================================================================
c=======================================================================
c     BEGIN: LOGML SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TIME SERIE (Y[T]:T=1,...,N), A SET OF CLUSTER LABELS
c     (E[T]:T=1,...,N) WITH 1=E[1]<=E[2]<=···<=E[N], AND
c     J IN {1,...,E[N]}, "LOGML" RETURNS THE LOG MARGINAL LIKELIHOOD
c     ML(Y[T]:T IN S[J]|THETA) UNDER DIFFERENT STATISTICAL MODELS.
c     HERE, THETA IS A PARAMETER THAT CONTROLS ML(·) AND
c     S[J]={I IN {1,...,N}:E[I]=J}.
c=======================================================================
c     INPUTS
c=======================================================================
c     indtype: TYPE OF MARGINAL LIKELIHOOD TO BE EVALUATED
      integer indtype
c     nobs: LENGTH OF THE TIME SERIES (Y[T]:T=1,...,N)
      integer nobs
c     obs: TIME SERIES (Y[T]:T=1,...,N)
      real*8 obs(nobs)
c     npars: MAXIMUM NUMBER OF ADMISSIBLE PARAMETERS
      integer npars
c     pars: PARAMETER THETA
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
c     OTHERS
      real*8 rw
c=======================================================================
c     ALGORITHM
c=======================================================================
      if (indtype.eq.1) then
         call lognornig(nobs,obs,npars,pars,labels,indj,rw)
         val=rw
      else if (indtype.eq.2) then
         call logpoigam(nobs,obs,npars,pars,labels,indj,rw)
         val=rw
      else
         print *,'Error in logml subroutine: check indtype'
         stop
      end if
c=======================================================================
      return
c     END: LOGML SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine rwmhp(noldps,oldps,vars,bins,nu,mu,invsigma,newps)
c=======================================================================
c=======================================================================
c     BEGIN: RWMHP SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO VECTORS (P^{OLD}[I]:I=1,...,M) IN (0,1)^M AND
c     (B[I]:I=1,...,M) IN {0,1}^M, "RWMHP" SEQUENTIALLY GENERATES
c     PROPOSALS P^{NEW}[I]~NORMAL(P^{OLD}[I],V[I]) AND ACCEPT THEM
c     (SET P^{OLD}[I]=P^{NEW}[I]) WITH PROBABILITY
c
c     MIN{1,F[I](P^{NEW}[I]|P^{OLD}[-I])/F[I](P^{OLD}[I]|P^{OLD}[-I])}.
c
c     NOTICE THAT P^{NEW}[I] IS IMMEDIATELY REJECTED IF P^{NEW}[I]
c     NOT IN (0,1). ON THE OTHER HAND, EACH FULL CONDITIONAL DENSITY
c     F[I](P|P^{OLD}[-I]) WITH P IN (0,1) IS PROPORTIONAL TO
c
c     (P^B[I])*{(1-P)^(1-B[I])}*LOGIT-T[M](P*|NU,MU,SIGMA),
c
c     WHERE
c
c     P*=(P^{OLD}[1],...,P^{OLD}[I-1],P,P^{OLD}[I+1],...,P^{OLD}[M])
c
c     AND LOGIT-T[M](·|NU,MU,SIGMA) IS THE M-DIMENSIONAL LOGIT-T
c     DENSITY WITH DEGREES OF FREEDOM NU, LOCATION VECTOR MU AND
c     SCALE MATRIX SIGMA.
c=======================================================================
c     INPUTS
c=======================================================================
c     noldps: LENGTH OF VECTOR (P^{OLD}[I]:I=1,...,M)
      integer noldps
c     pold: VECTOR (P^{OLD}[I]:I=1,...,M)
      real*8 oldps(noldps)
c     vars: VARIANCES FOR SYMMETRIC PROPOSALS (V[I]:I=1,...,M)
      real*8 vars(noldps)
c     bins: BINARY VECTOR (B[I]:I=1,...,M)
      integer bins(noldps)
c     nu: DEGREES OF FREEDOM NU
      real*8 nu
c     mu: LOCATION VECTOR MU
      real*8 mu(noldps)
c     sigma: INVERSE OF SCALE MATRIX SIGMA
      real*8 invsigma(noldps,noldps)
c=======================================================================
c     OUTPUTS
c=======================================================================
c     newps: PROBABILITY PARAMETERS (P^{NEW}[I]:I=1,...,M)
      real*8 newps(noldps)
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
      real*8 normr
      real*8 unifr
c=======================================================================
c     FORTRAN SUBROUTINES
c=======================================================================
c     mvtlogitd(···)
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
      integer jj
c     EXCLUSIVE FOR UPDATING (P^{OLD}[I]:I=1,...,M)
      real*8 old(noldps)
      real*8 new(noldps)
      real*8 logratio
c     OTHERS
      real*8 bw
      real*8 rw
      real*8 sw
      real*8 uw
c=======================================================================
c     ALGORITHM
c=======================================================================
      do ii=1,noldps
         old(ii)=oldps(ii)
         new(ii)=oldps(ii)
      end do
      do ii=1,noldps
c        P^{NEW}[I]~NORMAL(P^{OLD}[I],V[I])
         rw=normr(old(ii),vars(ii))
c        FIRST CONDITION FOR ACCEPTING P^{NEW}[I]:
c        · 0<P^{NEW}[I]<1
         if ((rw.gt.0.d0).and.(rw.lt.1.d0)) then
            new(ii)=rw
c           BW=B[I]
            bw=dble(bins(ii))
c           RW=LOGIT-T[M](P^{NEW}|NU,MU,SIGMA)
            call mvtlogitd(noldps,new,nu,mu,invsigma,rw)
c           SW={B[I]*LOG(P^{NEW}[I])}+{(1-B[I])*LOG(1-P^{NEW}[I])}
            sw=(bw*dlog(new(ii)))+((1.d0-bw)*dlog(1.d0-new(ii)))
            logratio=rw+sw
c           RW=LOGIT-T[M](P^{OLD}|NU,MU,SIGMA)
            call mvtlogitd(noldps,old,nu,mu,invsigma,rw)
c           SW={B[I]*LOG(P^{OLD}[I])}+{(1-B[I])*LOG(1-P^{OLD}[I])}
            sw=(bw*dlog(old(ii)))+((1.d0-bw)*dlog(1.d0-old(ii)))
            logratio=logratio-(rw+sw)
c           SECOND CONDITION FOR ACCEPTING P^{NEW}[I]:
c           · LOGRATIO>LOG(U), WHERE U~U(0,1)
            uw=unifr(0.d0,1.d0)
            rw=dlog(uw)
            if (logratio.gt.rw) then
c              UPDATING P^{OLD}[I]
               newps(ii)=new(ii)
               old(ii)=new(ii)
            else
               newps(ii)=old(ii)
               new(ii)=old(ii)
            end if
         end if
      end do
c=======================================================================
      return
c     END: RWMHP SUBROUTINE
      end
c=======================================================================
c=======================================================================

