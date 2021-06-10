c=======================================================================
c=======================================================================
c     TOOLS: (G)IBBS (S)AMPLING
c
c     SUBROUTINE NAMES:
c     · LOWERINDEX
c     · RASSOCGSAS
c     · RASSOCGUNI
c     · RWEIGHTBET
c     · UPPERINDEX
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
c      subroutine lowerindex(ntclust,tclust,indt,indit)
c=======================================================================
c=======================================================================
c     BEGIN: LOWERINDEX SUBROUTINE
c      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TEMPORAL CLUSTERING C=(C[1],...,C[N-1]) OF {1,...,N}
c     AND 0<T<N, "LOWERINDEX" RETURNS THE INDEX I[T] DEFINED AS
c
c     I[T]=
c
c     · MAX{S=1,...,T-1:C[S]=1}, IF T>1 AND
c       SUM(C[S]:S=1,...,T-1)>0.
c
c     · 0, OTHERWISE.
c=======================================================================
c     INPUTS
c=======================================================================
c     ntclust: LENGTH OF THE TEMPORAL AXIS (N)
c      integer ntclust
c     tclust: TEMPORAL CLUSTERING (C)
c      integer tclust(ntclust-1)
c     indt: POINT ON THE TEMPORAL AXIS (T)
c     integer indt
c=======================================================================
c     OUTPUTS
c=======================================================================
c     indit: INDEX (I[T])
c      integer indit
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
c     integer ss
c     OTHERS
c     integer rw
c=======================================================================
c     ALGORITHM
c=======================================================================
c     · CONDITION T=1 IMPLIES I[T]=0
c      if (indt.eq.1) then
c         indit=0
c      else
c        RW=SUM(C[S]:S=1,...,T-1)
c        rw=0
c         do ss=1,(indt-1)
c            rw=rw+tclust(ss)
c         end do
c        · CONDITION SUM(C[S]:S=1,...,T-1)=0 IMPLIES I[T]=0
c         if (rw.eq.0) then
c            indit=0
c         else
c           I[T]=MAX{S=1,...,T-1:C[S]=1}
c            indit=indt-1
c            ss=indt-1
c            do while (ss.gt.0)
c               if (tclust(ss).eq.0) then
c                  indit=indit-1
c                  ss=ss-1
c               else
c                  ss=0
c               end if
c            end do
c         end if
c      end if
c=======================================================================
c      return
c     END: LOWERINDEX SUBROUTINE
c      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine rassocgsas(ntclusts,tclusts,shpsa,shpsb,weight,
     & rassocg)
c=======================================================================
c=======================================================================
c     BEGIN: RASSOCGSAS SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO TEMPORAL CLUSTERINGS C[1],C[2] OF {1,...,N}, WHERE
c     C[L]=(C[1,L],...,C[N-1,L]), AND A VALUE W IN [0,1], "RASSOCGSAS"
c     RETURNS A RANDOM NUMBER G FROM PDF(G|C[1],C[2],W) INDUCED BY THE
c     FOLLOWING HIERARCHICAL MODEL:
c
c     A) PR(C[1],C[2]|G)=F(G)*PRODUCT(PR(C[L]):L=1,2), WHERE
c
c        F(G)=1+[G*PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]
c
c        AND
c
c        PR(C[L])=BE(A[L]+SIGMA[L],B[L]+(N-1)-SIGMA[L])/BE(A[L],B[L]).
c
c        · SHAPE PARAMETERS (A[L],B[L]:L=1,2) ARE FIXED.
c
c        · F(C[L])=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L]).
c
c        · SIGMA[L]=SUM(C[T,L]:T=1,...,N-1).
c
c        · BE(·,·): BETA FUNCTION.
c
c     B) G|W~{W*DELTA[0]}+{(1-W)*UNIFORM([-L[B],U[B]]-{0})}.
c         THIS IS A PRIOR WITH A DIRAC SPIKE AT 0 AND A UNIFORM SLAB
c         ON [-L[B],U[B]]-{0}.
c
c        · DELTA[0]: DIRAC MEASURE CENTRED AT 0.
c
c        · L[B]=F[B]/MAX{ALPHA[1]*ALPHA[2],BETA[1]*BETA[2]}].
c
c        · U[B]=F[B]/MAX{ALPHA[1]*BETA[2],ALPHA[2]*BETA[1]}].
c
c        · F[B]=(ALPHA[1]+BETA[1])*(ALPHA[2]+BETA[2]).
c
c     WITH THIS INFORMATION,
c
c     CDF(G|C[1],C[2],W)=
c
c       ·(W/[W+{(1-W)*F[1]}])*DELTA[0]((-INF,G])+
c         {(1-W)*F[1]/[W+{(1-W)*F[1]}]}*CDF(G|C[1],C[2]),
c         FOR ALL G IN (-INF,INF).
c
c     HERE, CDF(·|C[1],C[2]) IS THE CDF CORRESPONDING TO
c
c     PDF(G|C[1],C[2])=
c     
c       · {1+(G*F[2]})}/F[1], IF G IN [-L[B],U[B]]-{0}.
c
c       · 0, OTHERWISE.
c
c     0) F[1]=(U[B]+L[B])+{F[2]*[{(U[B]^2)-(L[B]^2)}/2]}.
c
c     1) F[2]=PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2).
c=======================================================================
c     INPUTS
c=======================================================================
c     ntclusts: LENGTH OF THE TEMPORAL AXIS (N)
      integer ntclusts
c     tclusts: TEMPORAL CLUSTERINGS (C[L]:L=1,2)
      integer tclusts(ntclusts-1,2)
c     shpsa: SHAPE PARAMETERS (A[L]:L=1,2)
      real*8 shpsa(2)
c     shpsb: SHAPE PARAMETERS (B[L]:L=1,2)
      real*8 shpsb(2)
c     weight: PROBABILITY WEIGHT (W)
      real*8 weight
c=======================================================================
c     OUTPUTS
c=======================================================================
c     rassocg: RANDOM NUMBER (G)
      real*8 rassocg
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
      real*8 unifr
c     - SEE FILE "TOOLSR.C" FOR MORE DETAILS
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ll,tt
c     EXCLUSIVE FOR STORING L[B]
      real*8 lowerb
c     EXCLUSIVE FOR STORING U[B]
      real*8 upperb
c     EXCLUSIVE FOR STORING SIGMA[L]
      real*8 sumc
c     OTHERS
      real*8 pw
      real*8 rw
      real*8 sw
      real*8 uw
c=======================================================================
c     ALGORITHM
c=======================================================================
      rw=(shpsa(1)+shpsb(1))*(shpsa(2)+shpsb(2))
      lowerb=rw/dmax1(shpsa(1)*shpsa(2),shpsb(1)*shpsb(2))
      upperb=rw/dmax1(shpsa(1)*shpsb(2),shpsa(2)*shpsb(1))
c     PW=PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]
c     · F(C[L])=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L])
      pw=1.d0
      do ll=1,2
         sumc=0.d0
         do tt=1,(ntclusts-1)
            sumc=sumc+dble(tclusts(tt,ll))
         end do
         rw=(shpsa(ll)+sumc)/(shpsa(ll)+dble(ntclusts-1)+shpsb(ll))
         pw=pw*(rw-(shpsa(ll)/(shpsa(ll)+shpsb(ll))))
      end do
c     SW=(W/[W+{(1-W)*F[1]}])
c     · F[1]=(U[B]+L[B])+{F[2]*[{(U[B]^2)-(L[B]^2)}/2]}
c     · F[2]=PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)
      rw=(upperb+lowerb)*(1.d0+(0.5d0*(upperb-lowerb)*pw))
      sw=weight/(weight+((1.d0-weight)*rw))
c     UW=UNIFORM(0,1)
      uw=unifr(0.d0,1.d0)
      if (uw.le.sw) then
         rassocg=0.d0
      else
         uw=unifr(0.d0,1.d0)
         if (pw.eq.0.d0) then
            rassocg=(-lowerb)+(uw*(lowerb+upperb))
         else
            rw=uw*((1.d0+(pw*upperb))**2)
            sw=(1.d0-uw)*((1.d0-(pw*lowerb))**2)
            rassocg=(-1.d0+dsqrt(rw+sw))/pw
         end if
      end if
c=======================================================================
      return
c     END: RASSOCGSAS SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine rassocguni(ntclusts,tclusts,shpsa,shpsb,rassocg)
c=======================================================================
c=======================================================================
c     BEGIN: RASSOCGUNI SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO TEMPORAL CLUSTERINGS C[1],C[2] OF {1,...,N}, WHERE
c     C[L]=(C[1,L],...,C[N-1,L]), "RASSOCGUNI" RETURNS A RANDOM
c     NUMBER G FROM PDF(G|C[1],C[2]) INDUCED BY THE FOLLOWING
c     HIERARCHICAL MODEL:
c
c     A) PR(C[1],C[2]|G)=F(G)*PRODUCT(PR(C[L]):L=1,2), WHERE
c
c        F(G)=1+[G*PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]
c
c        AND
c
c        PR(C[L])=BE(A[L]+SIGMA[L],B[L]+(N-1)-SIGMA[L])/BE(A[L],B[L]).
c
c        · SHAPE PARAMETERS (A[L],B[L]:L=1,2) ARE FIXED.
c
c        · F(C[L])=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L]).
c
c        · SIGMA[L]=SUM(C[T,L]:T=1,...,N-1).
c
c        · BE(·,·): BETA FUNCTION.
c
c     B) G~UNIFORM(-L[B],U[B]).
c
c        · L[B]=F[B]/MAX{ALPHA[1]*ALPHA[2],BETA[1]*BETA[2]}].
c
c        · U[B]=F[B]/MAX{ALPHA[1]*BETA[2],ALPHA[2]*BETA[1]}].
c
c        · F[B]=(ALPHA[1]+BETA[1])*(ALPHA[2]+BETA[2]).
c
c     WITH THIS INFORMATION,
c
c     PDF(G|C[1],C[2])=
c     
c       · {1+(G*F[1])}/F[2], IF G IN [-L[B],U[B]].
c
c       · 0, OTHERWISE.
c
c     0) F[1]=PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2).
c
c     1) F[2]=(U[B]+L[B])+{F[1]*[{(U[B]^2)-(L[B]^2)}/2]}.
c=======================================================================
c     INPUTS
c=======================================================================
c     ntclusts: LENGTH OF THE TEMPORAL AXIS (N)
      integer ntclusts
c     tclusts: TEMPORAL CLUSTERINGS (C[L]:L=1,2)
      integer tclusts(ntclusts-1,2)
c     shpsa: SHAPE PARAMETERS (A[L]:L=1,2)
      real*8 shpsa(2)
c     shpsb: SHAPE PARAMETERS (B[L]:L=1,2)
      real*8 shpsb(2)
c=======================================================================
c     OUTPUTS
c=======================================================================
c     rassocg: RANDOM NUMBER (G)
      real*8 rassocg
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
      real*8 unifr
c     - SEE FILE "TOOLSR.C" FOR MORE DETAILS
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ll,tt
c     EXCLUSIVE FOR STORING L[B]
      real*8 lowerb
c     EXCLUSIVE FOR STORING U[B]
      real*8 upperb
c     EXCLUSIVE FOR STORING SIGMA[L]
      real*8 sumc
c     OTHERS
      real*8 pw
      real*8 rw
      real*8 sw
      real*8 uw
c=======================================================================
c     ALGORITHM
c=======================================================================
      rw=(shpsa(1)+shpsb(1))*(shpsa(2)+shpsb(2))
      lowerb=rw/dmax1(shpsa(1)*shpsa(2),shpsb(1)*shpsb(2))
      upperb=rw/dmax1(shpsa(1)*shpsb(2),shpsa(2)*shpsb(1))
c     PW=PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]
c     · F(C[L])=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L])
      pw=1.d0
      do ll=1,2
         sumc=0.d0
         do tt=1,(ntclusts-1)
            sumc=sumc+dble(tclusts(tt,ll))
         end do
         rw=(shpsa(ll)+sumc)/(shpsa(ll)+dble(ntclusts-1)+shpsb(ll))
         pw=pw*(rw-(shpsa(ll)/(shpsa(ll)+shpsb(ll))))
      end do
c     UW=UNIFORM(0,1)
      uw=unifr(0.d0,1.d0)
      if (pw.eq.0.d0) then
         rassocg=(-lowerb)+(uw*(lowerb+upperb))
      else
         rw=uw*((1.d0+(pw*upperb))**2)
         sw=(1.d0-uw)*((1.d0-(pw*lowerb))**2)
         rassocg=(-1.d0+dsqrt(rw+sw))/pw
      end if
c=======================================================================
      return
c     END: RASSOCGUNI SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine rweightbet(assocg,shpa,shpb,rweight)
c=======================================================================
c=======================================================================
c     BEGIN: RWEIGHTBET SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A VALUE G IN [-L[B],U[B]], "RWEIGHTBET" RETURNS A RANDOM
c     NUMBER W FROM PDF(W|G) INDUCED BY THE FOLLOWING HIERARCHICAL
c     MODEL:
c
c     A) G|W~{W*DELTA[0]}+{(1-W)*UNIFORM([-L[B],U[B]]-{0})}.
c         THIS IS A PRIOR WITH A DIRAC SPIKE AT 0 AND A UNIFORM SLAB
c         ON [-L[B],U[B]]-{0}.
c
c        · DELTA[0]: DIRAC MEASURE CENTRED AT 0.
c
c        · L[B]=F[B]/MAX{ALPHA[1]*ALPHA[2],BETA[1]*BETA[2]}].
c
c        · U[B]=F[B]/MAX{ALPHA[1]*BETA[2],ALPHA[2]*BETA[1]}].
c
c        · F[B]=(ALPHA[1]+BETA[1])*(ALPHA[2]+BETA[2]).
c
c     B) W~BETA(A[W],B[W]).
c
c        · SHAPE PARAMETERS A[W] AND B[W] ARE FIXED.
c
c     WITH THIS INFORMATION,
c
c     PDF(W|G)=
c
c       · BETA(W|A[W]+1,B[W]), IF G=0.
c
c       · BETA(W|A[W],B[W]+1), IF G IN [-L[B],U[B]]-{0}.
c
c     HERE, BETA(·|A,B) IS THE BETA PDF WITH SHAPE PARAMETERS A AND B.
c=======================================================================
c     INPUTS
c=======================================================================
c     assocg: ASSOCIATION PARAMETER (G)
      real*8 assocg
c     shpa: SHAPE PARAMETER (A[W])
      real*8 shpa
c     shpb: SHAPE PARAMETER (B[W])
      real*8 shpb
c=======================================================================
c     OUTPUTS
c=======================================================================
c     rweight: RANDOM NUMBER (W)
      real*8 rweight
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
      real*8 betar
c     - SEE FILE "TOOLSR.C" FOR MORE DETAILS
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     OTHERS
      real*8 aw
      real*8 bw
c=======================================================================
c     ALGORITHM
c=======================================================================
      if (assocg.eq.0.d0) then
         aw=shpa+1.d0
         bw=shpb
         rweight=betar(aw,bw)
      else
         aw=shpa
         bw=shpb+1.d0
         rweight=betar(aw,bw)
      end if
c=======================================================================
      return
c     END: RWEIGHTBET SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
c      subroutine upperindex(ntclust,tclust,indt,indjt)
c=======================================================================
c=======================================================================
c     BEGIN: UPPERINDEX SUBROUTINE
c      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TEMPORAL CLUSTERING C=(C[1],...,C[N-1]) OF {1,...,N}
c     AND 0<T<N, "UPPERINDEX" RETURNS THE INDEX J[T] DEFINED AS
c
c     J[T]=
c
c     · MIN{S=T+1,...,N-1:C[S]=1}, IF T<N-1 AND
c       SUM(C[S]:S=T+1,...,N-1)>0.
c
c     · N, OTHERWISE.
c=======================================================================
c     INPUTS
c=======================================================================
c     ntclust: LENGTH OF THE TEMPORAL AXIS (N)
c      integer ntclust
c     tclust: TEMPORAL CLUSTERING (C)
c      integer tclust(ntclust-1)
c     indt: POINT ON THE TEMPORAL AXIS (T)
c      integer indt
c=======================================================================
c     OUTPUTS
c=======================================================================
c     indjt: INDEX (J[T])
c      integer indjt
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
c      integer ss
c     OTHERS
c      integer rw
c=======================================================================
c     ALGORITHM
c=======================================================================
c     · CONDITION T=N-1 IMPLIES J[T]=N
c      if (indt.eq.(ntclust-1)) then
c         indjt=ntclust
c      else
c        RW=SUM(C[S]:S=T+1,...,N-1)
c         rw=0
c         do ss=(indt+1),(ntclust-1)
c            rw=rw+tclust(ss)
c         end do
c        · CONDITION SUM(C[S]:S=T+1,...,N-1)=0 IMPLIES J[T]=N
c         if (rw.eq.0) then
c            indjt=ntclust
c         else
c           J[T]=MIN{S=T+1,...,N-1:C[S]=1}
c            indjt=indt+1
c            ss=indt+1
c            do while (ss.gt.0)
c               if (tclust(ss).eq.0) then
c                  indjt=indjt+1
c                  ss=ss+1
c               else
c                  ss=0
c               end if
c            end do
c         end if
c      end if
c=======================================================================
c      return
c     END: UPPERINDEX SUBROUTINE
c      end
c=======================================================================
c=======================================================================
