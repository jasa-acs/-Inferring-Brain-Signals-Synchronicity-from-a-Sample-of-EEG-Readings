#include "MIC.h"
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//
/***********************************************************************************/
/* vec2adj: vector to adjacency matrices                                           */
/***********************************************************************************/
//
arma::mat vec2adj(arma::vec const& cl){
  int n = cl.n_elem;
  arma::mat adj(n,n,fill::zeros);
  for(int n1=0; n1<n; n1++){
    for(int n2=n1+1; n2<n; n2++){
      if(cl(n1) == cl(n2)){
        adj(n1,n2) += 1;    adj(n2,n1) += 1;
      }
    }
  }
  return adj;
}
//
/***********************************************************************************/
/* C.I. functions on a posterior sample                                            */
/***********************************************************************************/
//
arma::colvec postCI(arma::colvec sample, double const& lvl){
  int const Q1 = (sample.n_elem*(1.0 - lvl)/2.0);
  int const Q2 = (sample.n_elem/2);
  int const Q3 = (sample.n_elem - Q1);

  arma::colvec qs(3);
  std::nth_element(sample.begin(),          sample.begin() + Q1, sample.end());
  std::nth_element(sample.begin() + Q1 + 1, sample.begin() + Q2, sample.end());
  std::nth_element(sample.begin() + Q2 + 1, sample.begin() + Q3, sample.end());
  qs(0) = sample(Q1-1);
  qs(1) = sample(Q2-1);
  qs(2) = sample(Q3-1);
  return qs;
}
//
/***********************************************************************************/
/* vec2ind: vector to indicator matrix                                             */
/***********************************************************************************/
//
arma::mat vec2ind(arma::vec const& S, int const &nK){
  int np = S.n_elem;
  arma::mat matS(nK,np,fill::zeros);
  for(int ic=0; ic<np; ic++){
    matS(S[ic]-1, ic) = (int) 1;
  }
  return matS;
}
