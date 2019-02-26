#include "MIC.h"
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//
/***********************************************************************************/
/* mrmultinom: matrix multinomial                                                  */
/***********************************************************************************/
// Tested
// arma::mat mrmultinom(arma::mat const &llm){
//   int ni = llm.n_cols;      int nc = llm.n_rows;
//   mat samp = arma::zeros(size(llm));
//   vec rvec = randu<vec>(ni);                 // Alternative approach
//   for(int jj=0; jj<ni; jj++){
//     colvec prob = arma::exp(llm.col(jj) - logsum(llm.col(jj)));
//     prob += 0.0001; prob /= arma::accu(prob); //num hack
//     // colvec prob = arma::exp(llm.col(jj));  // Less stable version, but faster.
//     // prob /= std::accumulate(prob.begin(), prob.end(), 0.0);
//     int ii = 0;
//     if(rvec[jj] == 0.0){
//       samp(0, jj) = (int) 1;
//       continue;
//     }
//     while(rvec[jj] > 0.0){
//       rvec[jj] -= prob[ii];
//       if (rvec[jj] <= 0.0) {
//         samp(ii,jj) = (int) 1;
//         break;
//       } else ii++;
//       if(ii >= nc){
//         samp(ii-1, jj) = (int) 1;
//         break;
//       }
//     }
//   }
//   return samp;
// }
//
//
/***********************************************************************************/
/* mrmultinom: matrix multinomial                                                  */
/***********************************************************************************/
// Tested
arma::mat mrmultinom0(arma::mat const &llm){
  int ni = llm.n_cols;      int nc = llm.n_rows;
  mat samp = arma::zeros(size(llm));
  vec rvec = randu<vec>(ni);                 // Alternative approach
  for(int jj=0; jj<ni; jj++){
    colvec prob = arma::exp(llm.col(jj) - logsum(llm.col(jj)));
    // prob += 0.0001; prob /= arma::accu(prob); //num hack
    // colvec prob = arma::exp(llm.col(jj));  // Less stable version, but faster.
    // prob /= std::accumulate(prob.begin(), prob.end(), 0.0);
    int ii = 0;
    if(rvec[jj] == 0.0){
      samp(0, jj) = (int) 1;
      continue;
    }
    while(rvec[jj] > 0.0){
      rvec[jj] -= prob[ii];
      if (rvec[jj] <= 0.0) {
        samp(ii,jj) = (int) 1;
        break;
      } else ii++;
      if(ii >= nc){
        samp(ii-1, jj) = (int) 1;
        break;
      }
    }
  }
  return samp;
}
//
/***********************************************************************************/
/* sampcoh: sample coherence (TBeta)                                               */
/***********************************************************************************/
//
double sampcoh(arma::mat const &L, arma::mat const &C,
               double const &a, double const &b){
  mat tot = L % C;
  int p = L.n_cols;   int k = L.n_rows;
  double totsum = std::accumulate(tot.begin(), tot.end(), 0.0);
  int cts = 0;        double samp = 0.0;
  double aa = a + totsum; double bb = b + p - totsum;
  while(cts < 10 && samp<(1.0/k)){
    samp = Rf_rbeta(aa,bb);
    cts++;
  }
  if(samp < 1.0/k) samp = 1.0/k;
  return samp;
}
//
/***********************************************************************************/
/* rDir: Dirichlet random sampler                                                  */
/***********************************************************************************/
//
arma::vec rDir(arma::colvec const &tot, double const &dir0){
  int nc   = tot.n_elem;
  vec samp(nc);
  for(int i=0; i<nc; i++) samp[i] = Rf_rgamma(dir0+tot[i], 1.0);
  samp /= arma::accu(samp);
  return samp;
}
//
