#ifndef __MIC_H__
#define __MIC_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <stdio.h>    //cout and file output
#include <fstream>    //ofstream
#include <iostream>   //cout
#include <sstream>    //istringstream
#include <string>     //string

using namespace arma;
using namespace Rcpp;

/************************************************************************************/
/* Data Structure                                                                   */
/************************************************************************************/
struct Fdat{                          // FIELD DATA STRUCTURE------------------------;
  field<cube>    dta;                 // field of cubes on epoch lvl;
  int            ns;                  // Number of Subjects;
  vec            ne;                  // Number of Epochs per sub;
  int            np;                  // Number of units to be clustered;
  vec            d;                   // Dimensionality per sub;
};

/************************************************************************************/
/* Parameters Structure                                                             */
/************************************************************************************/
struct NIG{                           // Normal-iGamma Posterior---------------------;
  mat means;                          // Posterior sampled means  d*K;
  mat dcovs;                         // Posterior sampled stds   d*K;
  mat Emeans;                         // Posterior expected means d*K;
  mat Edcovs;                        // Posterior expected stds  d*K;
};
struct PAR{                           // Parameters to be tracked and updated--------;
  vec               pi;               // Mixing Probability;
  vec               alpha;            // Subject coherence to population;
  field<vec>        beta;             // Epochs coherence to subjects;
  field<cube>       L;                // Epochs clusters;
  cube              C;                // Subject clusters;
  mat               S;                // Population clusters;
  // double            ll_c;             // Pr(Y|theta, L);
  // double            ll_i;             // Pr(Y|theta);
  // double            ll_ei;            // Pr(Y|E(theta|Y));
  double            ll_e;             // Pr(Y|post.theta, post.L)
};
/************************************************************************************/
/* Prior Structure                                                                  */
/************************************************************************************/
struct prior{                         // Priors -------------------------------------;
  vec           dir0;                 // Dirichlet prior;
  field<mat>    mu0;                  // NIG means, initialized as epoch mean
  field<mat>    dcov0;                // NIG sigs, initialized as sample std
  double        b0;                   // NIG igamma prior (rate), E fixed as samp cov
  double        a1;                   // TBeta prior shape1
  double        b1;                   // TBeta prior shape2
};

struct subprior{                      // epoch level NIG prior  ---------------------;
  vec           mu0;
  vec           dcov0;
  double        b0;
};
/************************************************************************************/
/* Random Number Generators                                                         */
/************************************************************************************/
arma::mat     mrmultinom(arma::mat const &llm);
arma::mat     mrmultinom0(arma::mat const &llm);
double        sampcoh(arma::mat const &L, arma::mat const &C,
                            double const &a, double const &b);
arma::vec     rDir(arma::colvec const &tot, double const &dir0);
/************************************************************************************/
/* Utility Functions                                                                */
/************************************************************************************/
arma::rowvec  mixpr(arma::vec const &p, double const &coh);
double        logsum(arma::colvec const &prcol);
void          clustalign(arma::mat &now, arma::mat const &ref);
NIG           NIGpost(arma::mat const &data, subprior const &pr, arma::mat const &C);
arma::mat     nu_est(arma::mat const &C, double const &a);

/************************************************************************************/
/* Posterior summary utility                                                        */
/************************************************************************************/
arma::mat     vec2adj(arma::vec const& cl);
arma::colvec  postCI(arma::colvec sample, double const& lvl);
arma::mat     vec2ind(arma::vec const& S, int const &nK);

#endif
