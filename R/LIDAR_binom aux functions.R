#this function generates vmat, which is then used to generate the theta matrix
get.theta=function(nlk,gamma,ncomm,nloc){
  vmat=matrix(NA,nloc,ncomm)
  for (i in 1:(ncomm-1)){
    if (i==(ncomm-1)) cumsoma=nlk[,ncomm]
    if (i< (ncomm-1)) cumsoma=rowSums(nlk[,(i+1):ncomm])
    vmat[,i]=rbeta(nloc,nlk[,i]+1,cumsoma+gamma)
  }
  vmat[,ncomm]=1
  convertVtoTheta(vmat,rep(1,nloc))
}
