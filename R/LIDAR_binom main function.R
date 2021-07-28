LIDAR_binom=function(y,n,ncomm,a.phi,b.phi,gamma,ngibbs,nburn,
                     theta.post,phi.post){
  #useful stuff
  nloc=nrow(y)
  nspp=ncol(y)
  hi=0.999999
  lo=0.000001
  NminusY=n-y
  
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(0.5,ncomm,nspp)
  
  array.LSKP=array(NA,dim=c(nloc,nspp,ncomm,2))
  prob1=rep(1/ncomm,ncomm)
  for (i in 1:nloc){
    for (j in 1:nspp){
      for (oo in 1:2){
        if (oo==1) n1=n[i,j]-y[i,j]
        if (oo==2) n1=y[i,j]
        array.LSKP[i,j,,oo]=rmultinom(1,size=n1,prob=prob1)
      }
    }
  }

  #to store outcomes from gibbs sampler
  if (theta.post)  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  if (!theta.post) theta.out=matrix(0,nloc,ncomm)
  if (phi.post)  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  if (!phi.post) phi.out=matrix(0,ncomm,nspp)
  
  llk=rep(NA,ngibbs)
  
  #run gibbs sampler
  options(warn=2)
  zeroes=array(0,dim=c(nloc,nspp,ncomm))
  for (i in 1:ngibbs){
    print(i)   

    #sample z when y=0
    tmp0=samplez0(theta=theta, OneMinusPhi=1-phi, 
                  NminusY=NminusY, ncommun=ncomm, nloc=nloc, nspp=nspp,
                  zeroes=zeroes)
    array.LSKP[,,,1]=tmp0$ArrayLSK
    nlk0=tmp0$nlk
    nks0=tmp0$nks
    
    #sample z when y=1
    tmp1=samplez1(theta=theta, phi=phi, 
                  y=y, ncommun=ncomm, nloc=nloc, nspp=nspp,
                  zeroes=zeroes)
    array.LSKP[,,,2]=tmp1$ArrayLSK
    nlk1=tmp1$nlk
    nks1=tmp1$nks

    #get parameters  
    theta=get.theta(nlk=nlk0+nlk1,gamma,ncomm,nloc) #theta.true#
    theta[theta>hi]=hi; theta[theta<lo]=lo
    phi=matrix(rbeta(nspp*ncomm,nks1+a.phi,nks0+b.phi),ncomm,nspp) #phi.true#
    phi[phi>hi]=hi; phi[phi<lo]=lo
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo

    #re-order groups
    if (i%%50==0 & i<nburn){
      med=apply(theta,2,mean)
      ordem=order(med,decreasing=T)
      theta=theta[,ordem]
      phi=phi[ordem,]
      array.LSKP[,,,1]=array.LSKP[,,ordem,1]
      array.LSKP[,,,2]=array.LSKP[,,ordem,2]
    }
    
    #calculate logl and store results  
    llk[i]=sum(dbinom(y,size=n,prob=prob,log=T))
    
    if (theta.post) theta.out[i,]=theta
    if (!theta.post & i>=nburn) theta.out=theta.out+theta
    if (phi.post) phi.out[i,]=phi
    if (!phi.post & i>=nburn) phi.out=phi.out+phi
  }
  seq1=nburn:ngibbs
  nseq1=length(seq1)
  res=list(llk=llk)
  if (theta.post)  res$theta=theta.out[seq1,]
  if (!theta.post) res$theta=theta.out/nseq1
  if (phi.post)    res$phi=phi.out[seq1,]
  if (!phi.post)   res$phi=phi.out/nseq1
  
  res  
}