## usethis namespace: start
#' @useDynLib LidarLDA, .registration = TRUE
## usethis namespace: end
NULL

#' Fit LDA model to LIDAR data
#'
#' This function uses a Gibbs sampler to fit a LDA model to LIDAR data. These data
#' consist of two P x H matrices (matrices y and n), where P is for pixel and H is for
#' height class.
#'
#' @param y P x H matrix containing the number of returns for each pixel in each height class
#' @param n P x H matrix containing the number of incoming light pulses for each pixel in each height class
#' @param nclust maximum number of clusters
#' @param a.phi parameter > 0 for the prior of phi
#' @param b.phi parameter > 0 for the prior of phi
#' @param gamma parameter between 0 and 1 for the prior of the truncated stick-breaking prior
#' @param ngibbs number of iterations for the MCMC algorithm
#' @param nburn number of iterations to discard as burn in
#' @param theta.post should samples from the posterior distribution for theta be returned (TRUE or FALSE)? If FALSE, just the posterior mean is returned
#' @param phi.post should samples from the posterior distribution for phi be returned (TRUE or FALSE)? If FALSE, just the posterior mean is returned
#' @param theta.init initial values, if available, for the P x nclust theta matrix. Default value for theta.init is NULL.
#' @param phi.init initial values, if available, for the nclust x H phi matrix. Default value for phi.init is NULL.
#'
#' @return This function returns a list containing several elements:
#'               \itemize{
#'                  \item llk:  log-likelihood for each iteration. This is a vector of size ngibbs.
#'                  \item theta: estimated relative abundance of each cluster in each pixel.
#'                  If theta.post==T, then this consists of a matrix where rows are samples
#'                  from the posterior distribution.
#'                  If theta.post==F, then this consists of a P x K matrix (K is the number of
#'                  clusters and P is the number of pixels) containing the posterior mean.
#'                  \item phi: estimated absorptance probability in each height class for each cluster.
#'                  If phi.post==T, then this consists of a matrix where rows are samples
#'                  from the posterior distribution.
#'                  If phi.post==F, then this consists of a K x H matrix (K is the number
#'                  of clusters and H is the number of height classes) containing the posterior mean.
#'               }
#' @importFrom stats rbeta dbinom rmultinom
#' @export

LidarLDA=function(y,n,nclust,a.phi,b.phi,gamma,ngibbs,nburn,
                  theta.post,phi.post,
                  theta.init=NULL,phi.init=NULL){
  #useful stuff
  npix=nrow(y)
  nheight=ncol(y)
  hi=0.999999
  lo=0.000001
  NminusY=n-y

  #initial values
  if (is.null(theta.init))  theta=matrix(1/nclust,npix,nheight)
  if (!is.null(theta.init)) theta=theta.init
  if (is.null(phi.init))    phi=matrix(0.5,nclust,nheight)
  if (!is.null(phi.init))   phi=phi.init

  array.LSKP=array(NA,dim=c(npix,nheight,nclust,2))
  prob1=rep(1/nclust,nclust)
  for (i in 1:npix){
    for (j in 1:nheight){
      for (oo in 1:2){
        if (oo==1) n1=n[i,j]-y[i,j]
        if (oo==2) n1=y[i,j]
        array.LSKP[i,j,,oo]=stats::rmultinom(1,size=n1,prob=prob1)
      }
    }
  }

  #to store outcomes from gibbs sampler
  if (theta.post)  theta.out=matrix(NA,ngibbs,nclust*npix)
  if (!theta.post) theta.out=matrix(0,npix,nclust)
  if (phi.post)  phi.out=matrix(NA,ngibbs,nclust*nheight)
  if (!phi.post) phi.out=matrix(0,nclust,nheight)

  llk=rep(NA,ngibbs)

  #run gibbs sampler
  options(warn=2)
  zeroes=array(0,dim=c(npix,nheight,nclust))
  for (i in 1:ngibbs){
    print(i)

    #sample z when y=0
    tmp0=samplez0(theta=theta, OneMinusPhi=1-phi,
                  NminusY=NminusY, nclust=nclust, npix=npix, nheight=nheight,
                  zeroes=zeroes)
    array.LSKP[,,,1]=tmp0$ArrayLSK
    nlk0=tmp0$nlk
    nks0=tmp0$nks

    #sample z when y=1
    tmp1=samplez1(theta=theta, phi=phi,
                  y=y, nclust=nclust, npix=npix, nheight=nheight,
                  zeroes=zeroes)
    array.LSKP[,,,2]=tmp1$ArrayLSK
    nlk1=tmp1$nlk
    nks1=tmp1$nks

    #get parameters
    theta=get.theta(nlk=nlk0+nlk1,gamma,nclust,npix) #theta.true#
    theta[theta>hi]=hi; theta[theta<lo]=lo
    phi=matrix(stats::rbeta(nheight*nclust,nks1+a.phi,nks0+b.phi),nclust,nheight) #phi.true#
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
    llk[i]=sum(stats::dbinom(y,size=n,prob=prob,log=T))

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
