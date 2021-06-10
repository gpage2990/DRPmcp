## Fortran subroutines

## Compile: Fortran codes
## R CMD SHLIB File1.f  ... FileN.f  File1.c ... FileN.c
## -o FileName
## If ToolsLA.f is needed: add -llapack -lblas after FileName


## "MCMCLOGITP" subroutine
mcmc.logitp = function(nburn,nskip,nsave,ydata,nu0,mu0,sigma0,mltype,
                      theta,rwmhvars){

  ndata = dim(ydata)[2]
  nseries = dim(ydata)[1]
  invsigma0 = solve(sigma0)
  nthetas = 4
  thetas = thetas[,1:nthetas]
  mcmcc = rep(0,nseries*((ndata-1)*nsave))
  mcmcp = rep(0,nseries*((ndata-1)*nsave))

  foo = .Fortran("MCMCLOGITP",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 nseries=as.integer(nseries),
                 ydata=as.double(ydata),
                 nu0=as.double(nu0),
                 mu0=as.double(mu0),
                 invsigma0=as.double(invsigma0),
                 mltypes=as.integer(mltypes),
                 nthetas=as.integer(nthetas),
                 thetas=as.double(thetas),
                 devs=as.double(devs),
                 mcmcc=as.integer(mcmcc),
                 mcmcp=as.double(mcmcp))

  C = array(foo$mcmcc,dim=c(nseries,ndata-1,nsave))
  P = array(foo$mcmcp,dim=c(nseries,ndata-1,nsave))

  return(list(C=C,P=P))
}



## "NORNIGNORNIGSAS" subroutine
nornignornigsas = function(nburn,nskip,nsave,ydata,alphas,betas,
                           c0,d0,mu0s,k0s,a0s,b0s,independent=FALSE)
{
  ndata = dim(ydata)[1]
  mctclusters = rep(0,2*(ndata-1)*nsave)
  mcnblocks = rep(0,2*nsave)
  mcassocg = rep(0,nsave)
  mcsasweight = rep(0,nsave)
  foo = .Fortran("NORNIGNORNIGSAS",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 ydata=as.double(ydata),
                 alphas=as.double(alphas),
                 betas=as.double(betas),
                 c0=as.double(c0),
                 d0=as.double(d0),
                 mu0s=as.double(mu0s),
                 k0s=as.double(k0s),
                 a0s=as.double(a0s),
                 b0s=as.double(b0s),
                 mctclusters=as.integer(mctclusters),
                 mcnblocks=as.integer(mcnblocks),
                 mcassocg=as.double(mcassocg),
                 mcsasweight=as.double(mcsasweight),
                 independent=as.integer(independent)
  )
  tclusters = matrix(foo$mctclusters,nrow=nsave,
                     ncol=(2*(ndata-1)),byrow=TRUE)
  nblocks = matrix(foo$mcnblocks,nrow=nsave,ncol=2,byrow=TRUE)
  assocg = c(foo$mcassocg)
  sasweight = c(foo$mcsasweight)
  return(list(tclusters=tclusters,nblocks=nblocks,assocg=assocg,
              sasweight=sasweight))
}


## "NORNIGNORNIGUNI" subroutine
nornignorniguni = function(nburn,nskip,nsave,ydata,alphas,betas,
                           mu0s,k0s,a0s,b0s)
{
  ndata = dim(ydata)[1]
  mctclusters = rep(0,2*(ndata-1)*nsave)
  mcnblocks = rep(0,2*nsave)
  mcassocg = rep(0,nsave)
  foo = .Fortran("NORNIGNORNIGUNI",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 ydata=as.double(ydata),
                 alphas=as.double(alphas),
                 betas=as.double(betas),
                 mu0s=as.double(mu0s),
                 k0s=as.double(k0s),
                 a0s=as.double(a0s),
                 b0s=as.double(b0s),
                 mctclusters=as.integer(mctclusters),
                 mcnblocks=as.integer(mcnblocks),
                 mcassocg=as.double(mcassocg)
                 )
  tclusters = matrix(foo$mctclusters,nrow=nsave,
                     ncol=(2*(ndata-1)),byrow=TRUE)
  nblocks = matrix(foo$mcnblocks,nrow=nsave,ncol=2,byrow=TRUE)
  assocg = c(foo$mcassocg)
  return(list(tclusters=tclusters,nblocks=nblocks,assocg=assocg))
}


## "NORNIGPOIGAMSAS" subroutine
nornigpoigamsas = function(nburn,nskip,nsave,ydata,alphas,betas,
                           c0,d0,mu0,k0,a0s,b0s)
{
  ndata = dim(ydata)[1]
  mctclusters = rep(0,2*(ndata-1)*nsave)
  mcnblocks = rep(0,2*nsave)
  mcassocg = rep(0,nsave)
  mcsasweight = rep(0,nsave)
  foo = .Fortran("NORNIGPOIGAMSAS",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 ydata=as.double(ydata),
                 alphas=as.double(alphas),
                 betas=as.double(betas),
                 c0=as.double(c0),
                 d0=as.double(d0),
                 mu0=as.double(mu0),
                 k0=as.double(k0),
                 a0s=as.double(a0s),
                 b0s=as.double(b0s),
                 mctclusters=as.integer(mctclusters),
                 mcnblocks=as.integer(mcnblocks),
                 mcassocg=as.double(mcassocg),
                 mcsasweight=as.double(mcsasweight)
  )
  tclusters = matrix(foo$mctclusters,nrow=nsave,
                     ncol=(2*(ndata-1)),byrow=TRUE)
  nblocks = matrix(foo$mcnblocks,nrow=nsave,ncol=2,byrow=TRUE)
  assocg = c(foo$mcassocg)
  sasweight = c(foo$mcsasweight)
  return(list(tclusters=tclusters,nblocks=nblocks,assocg=assocg,
              sasweight=sasweight))
}


## "NORNIGPOIGAMUNI" subroutine
nornigpoigamuni = function(nburn,nskip,nsave,ydata,alphas,betas,
                           mu0,k0,a0s,b0s)
{
  ndata = dim(ydata)[1]
  mctclusters = rep(0,2*(ndata-1)*nsave)
  mcnblocks = rep(0,2*nsave)
  mcassocg = rep(0,nsave)
  foo = .Fortran("NORNIGPOIGAMUNI",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 ydata=as.double(ydata),
                 alphas=as.double(alphas),
                 betas=as.double(betas),
                 mu0=as.double(mu0),
                 k0=as.double(k0),
                 a0s=as.double(a0s),
                 b0s=as.double(b0s),
                 mctclusters=as.integer(mctclusters),
                 mcnblocks=as.integer(mcnblocks),
                 mcassocg=as.double(mcassocg)
  )
  tclusters = matrix(foo$mctclusters,nrow=nsave,
                     ncol=(2*(ndata-1)),byrow=TRUE)
  nblocks = matrix(foo$mcnblocks,nrow=nsave,ncol=2,byrow=TRUE)
  assocg = c(foo$mcassocg)
  return(list(tclusters=tclusters,nblocks=nblocks,assocg=assocg))
}
