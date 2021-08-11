#this function generates vmat, which is then used to generate the theta matrix
get.theta=function(nlk,gamma,nclust,npix){
  vmat=matrix(NA,npix,nclust)
  for (i in 1:(nclust-1)){
    if (i==(nclust-1)) cumsoma=nlk[,nclust]
    if (i< (nclust-1)) cumsoma=rowSums(nlk[,(i+1):nclust])
    vmat[,i]=stats::rbeta(npix,nlk[,i]+1,cumsoma+gamma)
  }
  vmat[,nclust]=1
  convertVtoTheta(vmat,rep(1,npix))
}
