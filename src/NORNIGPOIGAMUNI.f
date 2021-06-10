c=======================================================================
c=======================================================================
      subroutine NORNIGPOIGAMUNI(nburn,nskip,nsave,ndata,ydata,
     & alphas,betas,mu0,k0,a0s,b0s,mctclusters,mcnblocks,mcassocg)
c=======================================================================
c=======================================================================
c     BEGIN: NORNIGPOIGAMUNI SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "NORNIGPOIGAMUNI" RETURNS POSTERIOR SAMPLES FROM THE FOLLOWING
c     BAYESIAN MODEL FOR CHANGE-POINT DETECTION:
c
c     A) PRIOR SPECIFICATION: TWO TEMPORAL CLUSTERINGS C[1],C[2] OF
c        {1,...,N}, WHERE C[L]=(C[1,L],...,C[N-1,L]) IN {0,1}^(N-1).
c
c     A.1) PR(C[1],C[2]|P[1],P[2])=
c          PRODUCT(BERN(C[T,L]|P[L]):T=1,...,N-1 AND L=1,2).
c
c          · BERN(·|P[L]): BERNOULLI PMF WITH PARAMETER P[L].
c
c     A.2) (P[1],P[2])|G FOLLOWS A CONTINUOUS DISTRIBUTION SUPPORTED
c          ON [0,1]^2, WITH PDF GIVEN BY
c
c          PDF(P[1],P[2]|G)=
c          W[G](P[1],P[2])*PRODUCT(BET(P[L]|ALPHA[L],BETA[L]):L=1,2),
c
c          · W[G](P[1],P[2])=
c            1+{G*PRODUCT(P[L]-{ALPHA[L]/(ALPHA[L]+BETA[L])}:L=1,2)}.
c
c          · G IN [-L[B],U[B]]: CONTROLS THE CORRELATION BETWEEN P[1]
c            AND P[2].
c
c          · BET(·|ALPHA[L],BETA[L]): BETA PDF WITH SHAPE PARAMETERS
c            ALPHA[L] AND BETA[L].
c
c          · SHAPES (A[L],B[L]:L=1,2) ARE FIXED.
c
c     A.3) G~UNIFORM(-L[B],U[B]).
c
c          · L[B]=F[B]*[MAX{ALPHA[1]*ALPHA[2],BETA[1]*BETA[2]}]^(-1).
c
c          · U[B]=F[B]*[MAX{ALPHA[1]*BETA[2],ALPHA[2]*BETA[1]}]^(-1).
c
c          · F[B]=(ALPHA[1]+BETA[1])*(ALPHA[2]+BETA[2]).
c
c     - SEE "LOGPR2YCF" SUBROUTINE IN FILE "TOOLSPM.F" FOR MORE
c       DETAILS ABOUT PR(C[1],C[2]|G).
c
c     B) DATA GENERATING MECHANISM: RECALL THAT EACH C[L] INDUCES A
c        UNIQUE INCREASING SEQUENCE OF B[L]+1 INTEGERS
c        (T[0,L],...,T[B[L],L]), WHERE B[L]=1+SUM(C[T,L]:T=1,...,N-1),
c        T[0,L]=0 AND
c
c        T[J,L]=
c
c          · MIN{S=T[J-1,L]+1,...,N-1:C[S,L]=1}, IF T[J-1,L]<N-1 AND
c            SUM(C[S,L]:S=T[J-1,L]+1,...,N-1)>0.
c
c          · N, OTHERWISE.
c
c        WITH THIS INFORMATION, C[L] REPRESENTS A UNIQUE CONTIGUOUS
c        SET PARTITION PI(C[L]) OF {1,...,N}, NAMELY,
c
c        PI(C[L])=CUP({T[J-1,L]+1,...,T[J,L]}:J=1,...,B[L]).
c
c     B.1) (Y[1,L],...,Y[N,L])|PI(C[L])~
c          PRODUCT(DF(Y[S,L]:S=T[J-1,L]+1,...,T[J,L]):J=1,...,B[L]).
c
c          · DF(·): DATA FACTOR.
c
c     C) DATA FACTOR MODEL: FOR ALL S IN {T[J-1,L]+1,...,T[J,L]} AND
c        J IN {1,...,B[L]},
c
c        Y[S,1]|MU[J,1],SIGMA[J,1]^2~NORMAL(MU[J,1],SIGMA[J,1]^2)
c
c        Y[S,2]|LAMBDA[J,2]~POISSON(LAMBDA[J,2]).
c
c        HERE, OBSERVATIONS ARE INDEPENDENT GIVEN MU[J,1],
c        SIGMA[J,1]^2 AND LAMBDA[J,2].
c
c     C.1) MU[J,1]|SIGMA[J,1]^2~NORMAL(MU[0,1],(SIGMA[J,1]^2)/K[0,1]).
c
c          · MEAN MU[0,1] AND PRECISION K[0,1] ARE FIXED.
c
c     C.2) SIGMA[J,1]^2~INVERSE-GAMMA(A[0,1],B[0,1]).
c
c          · SHAPE A[0,1] AND SCALE B[0,1] ARE FIXED.
c
c     C.3) LAMBDA[J,2]~GAMMA(A[0,2],B[0,2]).
c
c          · SHAPE A[0,2] AND RATE B[0,2] ARE FIXED.
c
c     - SEE "LOGDFNORNIG" AND "LOGDFPOIGAM" SUBROUTINES IN FILE
c       "TOOLSDF.F" FOR MORE DETAILS ABOUT
c       DF(Y[S,L]:S=T[J-1,L]+1,...,T[J,L]).
c=======================================================================
c     INPUTS: MCMC
c=======================================================================
c     nburn: BURN-IN ITERATIONS
      integer nburn
c     nskip: SKIP BETWEEN SAVED ITERATIONS
      integer nskip
c     nsave: SAVED ITERATIONS
      integer nsave
c=======================================================================
c     INPUTS: DATA
c=======================================================================
c     ndata: LENGTH OF THE TEMPORAL AXIS (N)
      integer ndata
c     ydata: TIME SERIES (Y[1],...,Y[N])
      real*8 ydata(ndata,2)
c=======================================================================
c     INPUTS: HYPER-PARAMETERS
c=======================================================================
c     alphas: SHAPES (ALPHA[L]:L=1,2)
      real*8 alphas(2)
c     betas: SHAPES (BETA[L]:L=1,2)
      real*8 betas(2)
c     mu0: MEAN (MU[0,1])
      real*8 mu0
c     k0: PRECISION (K[0,1])
      real*8 k0
c     a0s: SHAPES (A[0,L]:L=1,2)
      real*8 a0s(2)
c     b0s: SCALE AND RATE (B[0,L]:L=1,2)
      real*8 b0s(2)
c=======================================================================
c     OUTPUTS: MCMC SAMPLES
c=======================================================================
c     mctclusters: TEMPORAL CLUSTERINGS
      integer mctclusters(2*(ndata-1)*nsave)
c     mcnblocks: NUMBER OF BLOCKS
      integer mcnblocks(2*nsave)
c     mcassocg: ASSOCIATION PARAMETER
      real*8 mcassocg(nsave)
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
      real*8 unifr
c     - SEE FILE "TOOLSR.C" FOR MORE DETAILS
c=======================================================================
c     FORTRAN SUBROUTINES
c=======================================================================
c     logdfnornig(...)
c     - SEE FILE "TOOLSDF.F" FOR MORE DETAILS
c     logdfpoigam(...)
c     - SEE FILE "TOOLSDF.F" FOR MORE DETAILS
c     lowerindex(...)
c     - SEE FILE "TOOLSGS.F" FOR MORE DETAILS
c     rassocguni(...)
c     - SEE FILE "TOOLSGS.F" FOR MORE DETAILS
c     upperindex(...)
c     - SEE FILE "TOOLSGS.F" FOR MORE DETAILS
c     logpr2ycf(...)
c     - SEE FILE "TOOLSPM.F" FOR MORE DETAILS
c=======================================================================
c     WORKING VARIABLES: MCMC
c=======================================================================
c     actual: COUNTER (ONGOING ITERATION)
      integer actual
c     stored: COUNTER (STORED ITERATION)
      integer stored
c=======================================================================
c     WORKING VARIABLES: TEMPORAL CLUSTERING AND ASSOCIATION PARAMETER
c=======================================================================
c     tcluster: TEMPORAL CLUSTERINGS (C[L]:L=1,2)
      integer tclusters(ndata-1,2)
c     nblocks: NUMBER OF BLOCKS (B[L]:L=1,2)
      integer nblocks(2)
c     g: ASSOCIATION PARAMETER (G)
      real*8 g
c=======================================================================
c     WORKING VARIABLES: ALGORITHM
c=======================================================================
c     INDEXES
      integer l,s,t
      integer ll,tt
c     EXCLUSIVE FOR UPDATING C[1] AND C[2]
      integer tcluster(ndata-1)
      integer tcluster0(ndata-1,2)
      integer tcluster1(ndata-1,2)
      real*8 a0
      real*8 b0
      real*8 y0(ndata)
c     EXCLUSIVE FOR STORING I[T,L] AND J[T,L]
      integer it
      integer jt
c     EXCLUSIVE FOR STORING LOG RATIOS
      real*8 logratio1
      real*8 logratio2
      real*8 logvarpi
c     OTHERS
c      real*8 aw
c      real*8 bw
      real*8 rw
      real*8 uw
c=======================================================================
c     SETTING INITIAL VALUES
c=======================================================================
c     MCMC COUNTERS
      actual=1
      stored=1
c     NUMBER OF BLOCKS (B[L]=1:L=1,2)
      do l=1,2
         nblocks(l)=1
      end do
c     TEMPORAL CLUSTERINGS (C[T,L]=0:L=1,2 AND T=1,...,N-1)
      do l=1,2
         do t=1,(ndata-1)
            tclusters(t,l)=0
         end do
      end do
c     ASSOCIATION PARAMETER (G)
c      rw=(alphas(1)+betas(1))*(alphas(2)+betas(2))
c      aw=rw/dmax1(alphas(1)*alphas(2),betas(1)*betas(2))
c      bw=rw/dmax1(alphas(1)*betas(2),alphas(2)*betas(1))
c      uw=unifr(0.d0,1.d0)
c      g=(-aw)+(uw*(aw+bw))
      g=0.d0
c=======================================================================
c     METROPOLIS-HASTINGS-WITHIN-GIBBS ALGORITHM
c=======================================================================
c     PRINT ON SCREEN: BEGINING OF MCMC ITERATIONS
      print *,'Begining of MCMC iterations'
c=======================================================================
      call rndstart()
c     BEGIN: ITERATIONS
      do while (actual.le.(nburn+(nskip*nsave)))
c=======================================================================
c        UPDATING TEMPORAL CLUSTERINGS (C[L]:L=1,2)
c=======================================================================
         do l=1,2
c           A0=A[0,L]
            a0=a0s(l)
c           B0=B[0,L]
            b0=b0s(l)
c           Y0=(Y[T,L]:T=1,...,N)
            do t=1,ndata
               y0(t)=ydata(t,l)
            end do
c           UPDATING: C[T,L]:T=1,...,N-1
            do t=1,(ndata-1)
c              TCLUSTER0=(C[1],C[2]:C[T,L]=0)
c              TCLUSTER1=(C[1],C[2]:C[T,L]=1)
               do ll=1,2
                  do tt=1,(ndata-1)
                  tcluster0(tt,ll)=tclusters(tt,ll)
                  tcluster1(tt,ll)=tclusters(tt,ll)
                  end do
               end do
               tcluster0(t,l)=0
               tcluster1(t,l)=1
c              LOGRATIO1=LOG{PR(C[1],C[2]:C[T,L]=1|G)}
c              -LOG{PR(C[1],C[2]:C[T,L]=0|G)}
               call logpr2ycf(ndata,tcluster1,g,alphas,betas,rw)
               logratio1=rw
               call logpr2ycf(ndata,tcluster0,g,alphas,betas,rw)
               logratio1=logratio1-rw
c              TCLUSTER=C[L]
               do tt=1,(ndata-1)
                  tcluster(tt)=tclusters(tt,l)
               end do
c              IT=I[T,L]
               call lowerindex(ndata,tcluster,t,it)
c              JT=J[T,L]
               call upperindex(ndata,tcluster,t,jt)
c              LOGRATIO2=LOG{DF(Y[S,L]:S=I[T,L]+1,...,T)}
c              +LOG{DF(Y[S,L]:S=T+1,...,J[T,L])}
c              -LOG{DF(Y[S,L]:S=I[T,L]+1,...,J[T,L])}
               if (l.eq.1) then
                  call logdfnornig(ndata,y0,it,t,mu0,k0,a0,b0,rw)
                  logratio2=rw
                  call logdfnornig(ndata,y0,t,jt,mu0,k0,a0,b0,rw)
                  logratio2=logratio2+rw
                  call logdfnornig(ndata,y0,it,jt,mu0,k0,a0,b0,rw)
                  logratio2=logratio2-rw
               else
                  call logdfpoigam(ndata,y0,it,t,a0,b0,rw)
                  logratio2=rw
                  call logdfpoigam(ndata,y0,t,jt,a0,b0,rw)
                  logratio2=logratio2+rw
                  call logdfpoigam(ndata,y0,it,jt,a0,b0,rw)
                  logratio2=logratio2-rw
               end if
c              LOGVARPI=LOG(VARPI[T,L])
c              · LOG(VARPI[T,L])=LOG{PR(C[1],C[2]:C[T,L]=1|···)}
c                -LOG{PR(C[1],C[2]:C[T,L]=0|···)}
               logvarpi=logratio1+logratio2
c              UW=UNIFORM(0,1)
               uw=unifr(0.d0,1.d0)
               rw=dlog(uw)-dlog(1.d0-uw)
c              CONDITION:
c              · C[T,L]=1, IF LOG(VARPI[T,L])>LOG{U/(1-U)}
               if (logvarpi.gt.rw) then
                  tclusters(t,l)=1
               else
                  tclusters(t,l)=0
               end if
            end do
         end do
c=======================================================================
c        UPDATING NUMBER OF BLOCKS (B[L]:L=1,2)
c=======================================================================
         do l=1,2
c           S=SUM(C[T,L]:T=1,...,N-1)
            s=0
            do t=1,(ndata-1)
               s=s+tclusters(t,l)
            end do
            nblocks(l)=1+s
         end do
c=======================================================================
c        UPDATING ASSOCIATION PARAMETER (G)
c=======================================================================
         call rassocguni(ndata,tclusters,alphas,betas,rw)
         g=rw
c=======================================================================
c        PRINT ON SCREEN: BURN-IN PHASE COMPLETED
         if (actual.eq.nburn) then
            print *,'Burn-in phase completed'
         end if
c=======================================================================
c        PROCESSING: SAVED ITERATIONS
c        - AFTER NBURN ITERATIONS, NSAVE ITERATIONS BETWEEN NSKIP
c          STEPS
c=======================================================================
         if ((actual.gt.nburn).and.(mod(actual-nburn,nskip).eq.0)) then
c=======================================================================
c           STORING: MCMC SAMPLES
c=======================================================================
c           TEMPORAL CLUSTERINGS
            do l=1,2
               do t=1,(ndata-1)
                  s=t+((ndata-1)*(l-1))+(2*(ndata-1)*(stored-1))
                  mctclusters(s)=tclusters(t,l)
               end do
            end do
c           NUMBER OF BLOCKS
            do l=1,2
               s=l+(2*(stored-1))
               mcnblocks(s)=nblocks(l)
            end do
c           ASSOCIATION PARAMETER
            mcassocg(stored)=g
c=======================================================================
c           PRINT ON SCREEN: STORED ITERATION (MULTIPLES OF 100)
            if (mod(stored,100).eq.0) then
               print *,stored,'stored MCMC iterations'
            end if
c=======================================================================
c           UPDATING: COUNTER (STORED ITERATION)
c=======================================================================
            stored=stored+1
         end if
c=======================================================================
c        UPDATING: COUNTER (ONGOING ITERATION)
c=======================================================================
         actual=actual+1
      end do
c     END: MCMC ITERATIONS
      call rndend()
c=======================================================================
c     PRINT ON SCREEN: END OF MCMC ITERATIONS
      print *,'End of MCMC iterations'
c=======================================================================
      return
c     END: NORNIGPOIGAMUNI SUBROUTINE
      end
c=======================================================================
c=======================================================================
