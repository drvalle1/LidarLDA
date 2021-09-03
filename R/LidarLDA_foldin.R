#' Title
#'
#' @param y P x H matrix containing the number of returns for each pixel in each height class
#' @param n P x H matrix containing the number of incoming light pulses for each pixel in each height class
#' @param nclust maximum number of clusters (K)
#' @param gamma parameter between 0 and 1 for the prior of the truncated stick-breaking prior
#' @param ngibbs number of iterations for the MCMC algorithm
#' @param nburn number of iterations to discard as burn in
#' @param theta.post should samples from the posterior distribution for theta be returned (TRUE or FALSE)? If FALSE, just the posterior mean is returned
#' @param phi.post logical value (T or F) indicating if samples from the posterior distribution for phi are being provided or not
#' @param phi.estim if phi.post is true, then the algorithm
#' expects that this will be a large matrix containing
#' samples from the posterior distribution for phi. Each row of
#' this matrix should contain a separate sample. If phi.post is false,
#'  then this function expects that this is a K x H matrix
#'  containing an estimate (e.g., the posterior mean) of phi.
#'
#'
#' @return This function returns a list containing several elements:
#'               \itemize{
#'                  \item llk:  log-likelihood for each iteration. This is a vector of size ngibbs.
#'                  \item theta: estimated relative abundance of each cluster in each pixel.
#'                  This consists of a matrix where rows are samples
#'                  from the posterior distribution.
#'               }
#' @importFrom stats rbeta dbinom rmultinom
#' @export

LidarLDA_foldin=function(y,n,nclust,gamma,ngibbs,nburn,
                         phi.post,phi.estim,theta.post){
  #useful stuff
  npix=nrow(y)
  nheight=ncol(y)
  hi=0.999999
  lo=0.000001
  NminusY=n-y

  #initial values
  theta=matrix(1/nclust,npix,nclust)
  if (!phi.post) phi=phi.estim
  if (phi.post){ #do we have full posterior distribution
    npost=nrow(phi.estim)
    ooo=sample(npost,size=1)
    phi=matrix(phi.estim[ooo,],nclust,nheight)
  }

  array.LSKP=array(NA,dim=c(npix,nheight,nclust,2))
  for (i in 1:npix){
    for (j in 1:nheight){
      for (oo in 1:2){
        if (oo==1) {                #=0 (not observed)
          n1=n[i,j]-y[i,j]
          tmp=theta[i,]*(1-phi[,j])
        }
        if (oo==2) {                #=1 (observed)
          n1=y[i,j]
          tmp=theta[i,]*phi[,j]
        }
        prob1=tmp/sum(tmp)
        array.LSKP[i,j,,oo]=rmultinom(1,size=n1,prob=prob1)
      }
    }
  }
  # k=apply(array.LSKP,c(1,3),sum)
  # k1=k/apply(k,1,sum)
  # plot(theta.init,k1)
  # plot(theta,k1)

  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,nclust*npix)
  llk=rep(NA,ngibbs)

  #run gibbs sampler
  options(warn=2)
  zeroes=array(0,dim=c(npix,nheight,nclust))
  for (i in 1:ngibbs){
    print(i)

    #sample z when y=0
    tmp0=samplez0(theta=theta, OneMinusPhi=1-phi,
                  NminusY=NminusY, nclust=nclust, npix=npix,
                  nheight=nheight,zeroes=zeroes)
    array.LSKP[,,,1]=tmp0$ArrayLSK
    nlk0=tmp0$nlk
    # k=apply(array.LSKP[,,,1],c(1,3),sum)
    # unique(nlk0-k)

    #sample z when y=1
    tmp1=samplez1(theta=theta, phi=phi,
                  y=y, nclust=nclust, npix=npix, nheight=nheight,
                  zeroes=zeroes)
    array.LSKP[,,,2]=tmp1$ArrayLSK
    nlk1=tmp1$nlk
    # k=apply(array.LSKP[,,,2],c(1,3),sum)
    # unique(nlk1-k)

    #get parameters
    theta=get.theta(nlk=nlk0+nlk1,gamma,nclust,npix) #theta.true#
    # theta[theta>hi]=hi; theta[theta<lo]=lo

    if (phi.post){ #if we have full posterior distribution
      ooo=sample(npost,size=1)
      phi=matrix(phi.estim[ooo,],nclust,nheight)
    }

    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo

    #calculate logl and store results
    llk[i]=sum(dbinom(y,size=n,prob=prob,log=T))
    theta.out[i,]=theta
  }
  res=list(llk=llk)
  seq1=nburn:ngibbs
  if (theta.post)  res$theta=theta.out[seq1,]
  if (!theta.post) {
    tmp=colMeans(theta.out[seq1,])
    res$theta=matrix(tmp,npix,nclust)
  }
  res
}
