#include <RcppArmadillo.h>

using namespace Rcpp;

// This function makes multinomial draws
// [[Rcpp::export]]
IntegerVector rmultinom_1(NumericVector probs, int size) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(size, probs.begin(), k, ans.begin());
  return(ans);
}

// This function samples zs when ys=1
// [[Rcpp::export]]
List samplez1(NumericMatrix theta, NumericMatrix phi,
              IntegerMatrix y, int ncommun, int nloc, int nspp,
              NumericVector zeroes) {

  //convert array into arma::cube
  NumericVector vecArray=clone(zeroes);
  arma::cube ArrayLSK(vecArray.begin(),nloc, nspp, ncommun);

  IntegerMatrix nlk(nloc,ncommun);
  IntegerMatrix nks(ncommun,nspp);

  NumericVector prob(ncommun);
  IntegerVector znew(ncommun);

  for(int i=0; i<nloc; i++){
    for (int j=0; j<nspp; j++){
      for (int k=0; k<ncommun; k++){
        prob(k)=theta(i,k)*phi(k,j);
      }
      prob=prob/sum(prob);

      //multinomial draw
      znew=rmultinom_1(prob, y(i,j));

      for (int k=0; k<ncommun; k++){
        ArrayLSK(i,j,k)=znew[k];
      }

      //add to results to export
      nlk(i,_)=nlk(i,_)+znew;
      nks(_,j)=nks(_,j)+znew;
    }
  }
  return List::create(Named("nlk") = nlk,
                      Named("nks") = nks,
                      Named("ArrayLSK") = ArrayLSK);
}

// This function samples zs when ys=0
// [[Rcpp::export]]
List samplez0(NumericMatrix theta, NumericMatrix OneMinusPhi,
              IntegerMatrix NminusY, int ncommun, int nloc, int nspp,
              NumericVector zeroes) {

  //convert array into arma::cube
  NumericVector vecArray=clone(zeroes);
  arma::cube ArrayLSK(vecArray.begin(),nloc, nspp, ncommun);

  IntegerMatrix nlk(nloc,ncommun);
  IntegerMatrix nks(ncommun,nspp);

  NumericVector prob(ncommun);
  IntegerVector znew(ncommun);

  for(int i=0; i<nloc; i++){
    for (int j=0; j<nspp; j++){
      for (int k=0; k<ncommun; k++){
        prob(k)=theta(i,k)*OneMinusPhi(k,j);
      }
      prob=prob/sum(prob);

      //multinomial draw
      znew=rmultinom_1(prob, NminusY(i,j));

      for (int k=0; k<ncommun; k++){
        ArrayLSK(i,j,k)=znew[k];
      }

      //add to results to export
      nlk(i,_)=nlk(i,_)+znew;
      nks(_,j)=nks(_,j)+znew;
    }
  }
  return List::create(Named("nlk") = nlk,
                      Named("nks") = nks,
                      Named("ArrayLSK") = ArrayLSK);
}

// This function converts vmat into theta
// [[Rcpp::export]]
NumericMatrix convertVtoTheta(NumericMatrix vmat,
                              NumericVector prod) {
  NumericMatrix res(vmat.nrow(),vmat.ncol());

  for(int j=0; j<vmat.ncol();j++){
    res(_,j)=vmat(_,j)*prod;
    prod=prod*(1-vmat(_,j));
  }

  return (res);
}

// This function calculates ngreater
// [[Rcpp::export]]
IntegerMatrix ngreater(IntegerMatrix nlk,int nloc, int ncommun){
  IntegerMatrix ngreater(nloc,ncommun);
  int oo=ncommun-1;
  IntegerVector tmp(nloc);

  while (oo>=0){
    tmp=tmp+nlk(_,oo);
    ngreater(_,oo)=tmp;
    oo=oo-1;
  }
  return ngreater;
}
