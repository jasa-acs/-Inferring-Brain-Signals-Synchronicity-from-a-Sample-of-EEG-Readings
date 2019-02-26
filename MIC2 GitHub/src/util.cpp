#include "MIC.h"
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//
/***********************************************************************************/
/* mixpr: marginal mixing probability                                             */
/***********************************************************************************/
//
arma::rowvec mixpr(arma::vec const &p, double const &coh){
  int k = p.n_elem;
  double const tempc = (1.0 - coh) / (k - 1.0);
  arma::rowvec out(k);      out.fill(0.0);
  for(int kk = 0; kk<k; kk++) out(kk) = p(kk) * coh + (1.0-p(kk)) * tempc;
  return out;
}
//
/***********************************************************************************/
/* logsum                                                                          */
/***********************************************************************************/
//
double logsum(arma::colvec const &prcol){
  double m = prcol.min();       //int nl = prcol.n_elem;
  // vec temp(nl);
  // for(int ii=0; ii<nl; ii++) temp(ii) = std::exp(prcol(ii) - m);
  long double sum = arma::accu(arma::exp(prcol - m));
  // long double sum = std::accumulate(temp.begin(), temp.end(), 0.0);
  // double out = m + std::log (sum);
  // long double sum = arma::accu(temp);
  m += std::log(sum);
  return m;
}
//
/***********************************************************************************/
/* clustalign                                                                      */
/***********************************************************************************/
//' Cpp function for aligning indicator matrices
//'
//' \code{clustalign} performs indicator matrices alignment to achieve maximal
//'   concordance of \code{now} to the reference \code{ref}, by swapping the rows
//'   of matrix \code{now}. Please refer to \code{\link{align}} for general vector
//'   and matrix alignment.
//'
//' @param now,  matrix with its rows to be realigned
//' @param ref,  matrix to be used as reference
//' @return no value returned, and matrix \code{now} is aligned on the spot
//'
//[[Rcpp::export]]
void clustalign(arma::mat &now, arma::mat const &ref){
  int k = ref.n_rows;
  double best_score = arma::accu(now % ref);
  for(int c1=0; c1<k-1; c1++){
    for(int c2=c1+1; c2<k; c2++){
      now.swap_rows(c1,c2);
      double now_score = arma::accu(now % ref);
      if(now_score < best_score){
        now.swap_rows(c1,c2);
      } else {
        best_score = now_score;
      }
    }
  }
}
//
/***********************************************************************************/
/* NIGpost: Normal-inverseGamma posterior                                          */
/***********************************************************************************/
//
NIG NIGpost(arma::mat const &data, subprior const &pr, arma::mat const &C){
  NIG post;
  int d = data.n_rows;  int p = data.n_cols;  int k = C.n_rows;
  colvec datstd = var(data, 0, 1);
  post.means = randn<mat> (d,k); post.dcovs = zeros<mat> (d,k);
  post.Emeans= zeros<mat> (d,k); post.Edcovs= zeros<mat> (d,k);
  colvec csize = arma::sum(C,1);
  for(int c=0; c<k; c++){
    if(csize(c) < 1.0){ // empty cluster
      post.means.col(c)   = pr.mu0;
      post.dcovs.col(c)   = datstd;
      post.Emeans.col(c)  = pr.mu0;
      post.Edcovs.col(c)  = datstd;
      continue;
    }
    for(int dim=0; dim<d; dim++){
      // double shape = csize(c)/2.0 + 1.0 + (double) pr.b0/ pr.dcov0(dim); //fixed rate
      // double b = pr.b0;                                                  //fixed rate
      double shape = 2.0 + pr.b0 + csize(c) / 2.0;    //fixed shape
      double b     = (1.0 + pr.b0) * pr.dcov0(dim);   //fixed shape
      double smean=0.0;     double ssqr=0.0;
      for(int ep=0; ep<p; ep++){
        if(C(c,ep) > .9){
          smean += data(dim,ep);
          ssqr  += std::pow(data(dim,ep), 2.0);
        }
      }
      double mupost = (pr.mu0(dim) + smean) / (csize(c) + 1.0);
      b +=  0.5 * (ssqr - csize(c) * std::pow(smean/(csize(c)+0.0), 2.0) +
        std::pow((smean/(csize(c)+0.0) - pr.mu0(dim)), 2.0) *
        csize(c) / (csize(c) + 1.0));
      double tau = Rf_rgamma(shape, 1.0/b);
      post.dcovs(dim, c) = std::pow(tau, -1.0);
      double sig = std::pow(tau*(csize(c) + 1.0), -0.5);
      post.means(dim,c)  = post.means(dim,c) * sig + mupost;
      post.Emeans(dim,c) = mupost;
      post.Edcovs(dim,c) = b/(shape+1.0);
    }
  }
  return post;
}
//
/***********************************************************************************/
/* nu_est: connection calculation                                                  */
/***********************************************************************************/
//
arma::mat nu_est(arma::mat const &C, double const &a){
  double const tempcc = std::log(1.0 - a + 0.00001) - std::log(C.n_rows - 1.0);
  arma::mat nu = C * (std::log(a) - tempcc) + tempcc;
  return nu;
}

