#include <RcppArmadillo.h>

//' Compute the RV coefficient using Rcpp
//' 
//' @param mDx A numeric matrix of pairwise distances.
//' @param mDy A second numeric matrix of pairwise distances.
//' @param mC See the equation 2.4 in Josse and Homes manual.
//' 
//' @keywords Internal
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double RVcoeff(arma::mat mDx, arma::mat mDy, arma::mat mC){
  
  arma::mat Wx = mC*arma::pow(mDx,2)*mC; 
  arma::mat Wy = mC*arma::pow(mDy,2)*mC;
  
  // <Wx,Wy>
  double num = arma::trace(Wx*Wy);
  
  // ||Wx|| = sqrt(trace(Wx*Wx))
  double deno1 = arma::trace(Wx*Wx);
  double deno2 = arma::trace(Wy*Wy);
  
  
  
  double deno = sqrt(deno1*deno2); 
  
  double RV = num/deno;
  return RV;

  
}

/*** R
sourceCpp("RVcoeff.cpp")
 */
