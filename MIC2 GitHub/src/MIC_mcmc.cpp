#include "MIC.h"
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//
using namespace std;
/***********************************************************************************/
/* FILES & CONSTANTS ***************************************************************/
/***********************************************************************************/
//
#define fout_pop "post/pop.txt"
#define fout_sub "post/sub.txt"
#define fout_pre "post/sub_"
#define post_CIs "post/CIs.txt"
#define post_dis "post/dis.txt"
//
/***********************************************************************************/
/* STRUCTURES & GLOBAL VARIABLES ***************************************************/
/***********************************************************************************/
//
Fdat   dta;                   // Data Structure
prior   pr;                   // Prior parameters
PAR    par;                   // Parameters
//
// Output files ------------------------------------------------------------------
ofstream out1;                // population cluster (S); pi; ll's(3);
ofstream out2;                // subject clusters(Ci); alpha;
//
/***********************************************************************************/
/* OutputSample()                                                                  */
/***********************************************************************************/
//
void OutputSample()
{
  int K = par.pi.n_elem;
  {
    //--------------------------------------------------------------------------------
    // Population output;
    //--------------------------------------------------------------------------------
    std::stringstream ss;                             // population clusters
    for(int s=0; s<dta.np; s++){
      int outk=0; while(par.S(outk, s) < 0.9) outk++;
      ss << (int) outk+1 << " ";
    }
    for(int r=0; r<K; r++){                           // mixing probs *3
      ss << par.pi(r) << " ";
    }
    // ss << par.ll_c << " "; ss << par.ll_i << " "; ss << par.ll_ei << "\n"; //ll*3
    ss << "\n";   out1 << ss.str();
  }
  {
    //--------------------------------------------------------------------------------
    // Individual output;
    //--------------------------------------------------------------------------------
    std::stringstream ss;
    for(int is=0; is<dta.ns; is++){
      for(int ip=0; ip<dta.np; ip++){                 // individual clusters
        int outk=0; while(par.C.slice(is)(outk, ip) < 0.9) outk++;
        ss << (int) outk+1 << " ";
      }
      ss << par.alpha(is) << "\n";                    // alpha(sub)
    }
    out2 << ss.str();
  }
}
//
/***********************************************************************************/
/* PostSummary()                                                                   */
/***********************************************************************************/
//
void PostSummary(int const &nsims)
{
  // Constants acquisition
  int k       = par.pi.n_elem;            //number of clusters;
  int tot_ne  = accu(dta.ne);             //total number of epochs;
  int tot_c   = 1+dta.ns+tot_ne;          //total number of cluster vectors;
  int tot_ci  = k + dta.ns + tot_ne;      //total number of C.I. parameters;
  double du   = pow(dta.np, 2.0);         //upper limit of distance;
  //--------------------------------------------------------------------------------//
  // Read in Post Samples
  //--------------------------------------------------------------------------------//
  // Output files
  arma::mat dismat = zeros(nsims, tot_c);
  arma::mat cimat  = zeros(3, tot_ci);
  //
  //-----------------------pop.txt------------------------//
  {
    arma::mat avg(dta.np,dta.np);     avg.fill(0.0);
    {// To get mean adjacency
      std::string line;   std::ifstream in1(fout_pop);
      while (std::getline(in1, line)){
        std::istringstream iss(line);
        arma::vec clust(dta.np);
        for(int ini=0; ini<dta.np; ini++) iss >> clust[ini];
        avg += vec2adj(clust) / (nsims+0.0);
      }
      in1.close();
    }
    {// Get the S.best, clust_dist, pi's, ll's
      std::string line;   std::ifstream in1(fout_pop);
      int counter = 0;    arma::mat postpars(nsims, k);
      double best_d = du;
      while (std::getline(in1, line)){
        std::istringstream iss(line);
        arma::vec clust(dta.np);
        for(int cli=0; cli<dta.np; cli++) iss >> clust[cli];
        for(int cj=0; cj<k;cj++) iss >> postpars(counter, cj);
        double dist = arma::norm((vec2adj(clust) - avg), "fro");
        dismat(counter,0) = dist;
        if(dist < best_d) {
          par.S = vec2ind(clust, k);    best_d = dist;
        }
        counter++;
      }
      in1.close();
      arma::rowvec means = arma::mean(postpars, 0);
      for(int kk=0; kk<k; kk++) par.pi(kk) = means(kk);
      // par.ll_c = means(k);  par.ll_i = means(k+1); par.ll_ei = means(k+2);
      for(int cl=0; cl<k; cl++) cimat.col(cl) =
        postCI(postpars.col(cl), .90);
    }
  }
  //-----------------------sub.txt------------------------//
  {
    arma::cube avg(dta.np,dta.np,dta.ns);     avg.fill(0.0);
    {// To get mean adjacency cube
      std::string line;   std::ifstream in1(fout_sub);
      int counter = 0;
      while (std::getline(in1, line)){
        if(counter==dta.ns) counter=0; //reset
        std::istringstream iss(line);
        arma::vec clust(dta.np);
        for(int ini=0; ini<dta.np; ini++) iss >> clust[ini];
        avg.slice(counter) += vec2adj(clust) / (nsims+0.0);
        counter++;
      }
      in1.close();
    }
    {// Get the Ci.best, clust_dist, alpha's
      std::string line;   std::ifstream in1(fout_sub);
      int ct1   = 0;  arma::mat postpars(nsims, dta.ns);
      int ct2   = 0;  arma::vec best_d(dta.ns);   best_d.fill(du);
      while (std::getline(in1, line)){
        std::istringstream iss(line);   arma::vec clust(dta.np);
        for(int cli=0; cli<dta.np; cli++) iss >> clust[cli]; string rest;
        iss >> postpars(ct2, ct1);
        double dist = arma::norm((vec2adj(clust) - avg.slice(ct1)), "fro");
        dismat(ct2,ct1+1) = dist;
        if(dist < best_d[ct1]){
          best_d(ct1) = dist;   par.C.slice(ct1) = vec2ind(clust, k);
        }
        ct1++;  if(ct1==dta.ns) {ct2++;}
        ct1 %= dta.ns;
      }
      in1.close();
      arma::rowvec means = arma::mean(postpars, 0);
      for(int kk=0; kk<dta.ns; kk++) par.alpha(kk) = means(kk);
      for(int cl=0; cl<dta.ns; cl++) cimat.col(cl+k) =
        postCI(postpars.col(cl), .90);
    }
  }
  //----------------------sub_i.txt-----------------------//
  for(int sub=0; sub<dta.ns; sub++){
    std::string fname = fout_pre + std::to_string(sub) + ".txt";
    {// for each subject
      arma::cube avg(dta.np,dta.np,dta.ne(sub)); avg.fill(0.0);
      {// To get mean adjacency cube
        std::string line;   std::ifstream in1(fname);
        int counter = 0;
        while (std::getline(in1, line)){
          if(counter==dta.ne(sub)) counter=0; //reset
          std::istringstream iss(line);
          arma::vec clust(dta.np);
          for(int ini=0; ini<dta.np; ini++) iss >> clust[ini];
          avg.slice(counter) += vec2adj(clust) / (nsims+0.0);
          counter++;
        }
        in1.close();
      }
      //  recording offsets;
      int offset = 1 + dta.ns;
      for(int ss=0; ss<sub; ss++) offset += (int) dta.ne(ss);
      {// Get the L.best, clust_dist, beta's
        std::string line;   std::ifstream in2(fname);
        int ct1 = 0;    arma::mat postpars(nsims, dta.ne(sub));
        int ct2 = 0;    arma::vec best_d(dta.ne(sub)+1);
        best_d.fill(pow(dta.np, 2.0));
        while (std::getline(in2, line)){
          std::istringstream iss(line);   arma::vec clust(dta.np);
          if(ct1==dta.ne(sub)) ct1=0;
          for(int cli=0; cli<dta.np; cli++) iss >> clust[cli];
          iss >> postpars(ct2, ct1);
          double dist = arma::norm((vec2adj(clust) - avg.slice(ct1)), "fro");
          dismat(ct2,ct1+offset) = dist;
          if(dist < best_d(ct1)){
            best_d(ct1) = dist; par.L(sub).slice(ct1) = vec2ind(clust, k);
          }
          ct1++; if(ct1==dta.ne(sub)) ct2++;
        }
        arma::rowvec means = arma::mean(postpars, 0);
        for(int kk=0; kk<dta.ne(sub); kk++) par.beta(sub)(kk) = means(kk);
        for(int cl=0; cl<dta.ne(sub); cl++) cimat.col(cl+k+offset-1) =
          postCI(postpars.col(cl), .90);
        in2.close();
      }
    }
  }
  //-------------------OUTPUT-----------------------------//
  { // OUTPUT CI
    ofstream output;  output.open(post_CIs);
    for(int ni=0; ni<3; ni++){
      for(int nj=0; nj<tot_ci; nj++){
        output << cimat(ni, nj) << " ";
      }
      output << "\n";
    }
    output.close();
  }
  { // OUTPUT dist
    ofstream output;  output.open(post_dis);
    for(int ni=0; ni<nsims; ni++){
      std::stringstream ss;
      for(int nj=0; nj<tot_c; nj++){
        ss << dismat(ni,nj) << " ";
      }
      ss << "\n";       output << ss.str();
    }
    output.close();
  }
}
//
/***********************************************************************************/
/*  MIC_mcmc()                                                                     */
/***********************************************************************************/
//' MIC sampling function in C++
//'
//' \code{MIC_mcmc} performs Gibbs sampling in C++ as an essential part of
//'   \code{\link{MIC}}
//'
//' @param data, R list of 3D array, see \code{\link{MIC}}
//' @param K,    integer as number of clusters
//' @param run,  integer as number of MCMC samples
//' @param thin, integer as thinning every few iterations
//'
// [[Rcpp::export]]
List MIC_mcmc(Rcpp::List const &data,       // Data as R-List of 3D array:d,p,ne;
                        int const &K,       // Predetermined #of clusters;
                      int const &run,       // Iterations;
                     int const &thin       // thinning;
                )
  {
  // RNG from R set.seed() call
  RNGScope scope;
  // Reading in data
  dta.ns        = data.length();            // Number of subjects;
  dta.dta       = field<cube>(dta.ns);      // subfield of cubes;
  dta.ne        = vec(dta.ns);              // Number of epochs;
  dta.d         = vec(dta.ns);              // Dim per subject;
  for(int sub = 0; sub<dta.ns; sub++) {     // Reading in cubes
    cube cubedat = data[sub];
    dta.dta(sub) = cubedat;
    dta.np       = cubedat.n_cols;
    dta.d(sub)   = cubedat.n_rows;
    dta.ne(sub)  = cubedat.n_slices;
  }
  // parameter initialization                   ***    Values   ***
  par.pi = vec(K);                              par.pi.fill(1.0/K);
  par.S = mat(K, dta.np);                       par.S.fill(0);
  par.C = cube(K, dta.np, dta.ns);              par.C.fill(0);
  par.alpha = vec(dta.ns);                      par.alpha.fill(1.0/K);
  par.beta = field<vec>(dta.ns);
  par.L    = field<cube>(dta.ns);
  for(int sub=0; sub<dta.ns; sub++){
    par.beta(sub) = vec(dta.ne(sub));           par.beta(sub).fill(1.0/K);
    par.L(sub) = cube(K, dta.np, dta.ne(sub));  par.L(sub).fill(0);
  }
  // prior declaration:
  pr.b0   = 0.0001;                                      //NIG(uninformative)
  pr.a1   = 1;                   pr.b1   = 1;           //Tbeta
  pr.dir0 = vec(K);              pr.dir0.fill(3.0);     //Dirichlet********************
  pr.mu0  = field< mat >(dta.ns);                       //NIG-means
  pr.dcov0= field< mat >(dta.ns);                       //NIG-vars
  for(int sub=0; sub<dta.ns; sub++){
    pr.mu0(sub)   = mat(dta.d(sub),dta.ne(sub));      pr.mu0(sub).fill(0.0);
    pr.dcov0(sub) = mat(dta.d(sub),dta.ne(sub));      pr.dcov0(sub).fill(0.0);
  }
  // Initial Clusters, NIG priors and reference label:
  {
    urowvec init(dta.np);
    for(int sub=0; sub<dta.ns; sub++){
      for(int ep=0; ep<dta.ne(sub); ep++){
        mat matdata = dta.dta(sub).slice(ep);
        gmm_diag model;
        bool status = model.learn(matdata, K, maha_dist,
                                  random_spread, 10, 10, 1e-10, false);
        if(status != false){
          init = model.assign(matdata, prob_dist);
        }
        for(int unit=0; unit<dta.np; unit++){
          int initc = init(unit);
          par.L(sub)(initc, unit, ep) = 1;
        }
        // align clusters to the first epoch(slice)
        clustalign(par.L(sub).slice(ep), par.L(sub).slice(0));
        // NIG priors
        pr.mu0(sub).col(ep)   += arma::mean(model.means, 1);
        pr.dcov0(sub).col(ep) += arma::mean(model.dcovs, 1);
      }
    }
  }
  // Open output files ----------------------------------------------
  ofstream *out3 = new ofstream[dta.ns];              // sub_i.txt epoch output
  for(int s=0; s<dta.ns; s++) out3[s].open(fout_pre +
      std::to_string(s) + ".txt");
  out1.open(fout_pop);                                // pop.txt   pop output
  out2.open(fout_sub);                                // sub.txt   sub output
  // Output counter
  int oc = 0;
  // ----------------------------Gibbs Sampler---------------------------------------
  for(int iter=0; iter<run; iter++){
    par.pi /= std::accumulate(par.pi.begin(), par.pi.end(), 0.0);
    // double ll_c=0.0;    double ll_i=0.0;    double ll_ei=0.0; //ll*3
    // --------------------------Subject level---------------------------------------
    for(int sub=0; sub<dta.ns; sub++){
      rowvec mpi = mixpr(par.pi, par.alpha(sub));
      // --------------------------Epoch level---------------------------------------
      for(int ie=0; ie<dta.ne(sub); ie++){
        subprior spr;
        spr.mu0     = pr.mu0(sub).col(ie);      spr.dcov0  = pr.dcov0(sub).col(ie);
        spr.b0      = pr.b0;
        mat epdta = dta.dta(sub).slice(ie);
        // --------------------------------------------------------------------------
        // EPmodule 1: mixing prob | alpha, beta, Pi
        // --------------------------------------------------------------------------
        vec tempi = arma::conv_to< vec >::from(mpi);
        rowvec mpij = mixpr(tempi, par.beta(sub)(ie));
        mpij /= arma::accu(mpij);
        mpij[0] = 1.0 - arma::accu(mpij(span(1,K-1)));
        // --------------------------------------------------------------------------
        // EPmodule 2: NIG | L, data, pr
        // --------------------------------------------------------------------------
        NIG post = NIGpost(epdta, spr, par.L(sub).slice(ie));
        // --------------------------------------------------------------------------
        // EPmodule 3: gmm probs | post.mean, post.sigs, mpij
        // --------------------------------------------------------------------------
        arma::gmm_diag epmodel;
        epmodel.set_params(post.means, post.dcovs, mpij);
        mat ll = zeros(K, dta.np);
        for(int ek=0; ek<K; ek++){
          ll.row(ek) = epmodel.log_p(epdta, ek);
        }
        // mat ll2 = ll; //for IC calculation
        // if(iter > 0) ll += nu_est(par.C.slice(sub), par.beta(sub)(ie)); //MIC1
        if(iter > 0) ll.each_col() += arma::log(mpij.t());              //MIC2
        // --------------------------------------------------------------------------
        // EPmodule 4: L | gmm_probs
        // --------------------------------------------------------------------------
        par.L(sub).slice(ie) = mrmultinom0(ll);
        // Realign the clusters
        if (iter > 0) clustalign(par.L(sub).slice(ie), par.S);
        // --------------------------------------------------------------------------
        // EPmodule 5: ICs: ll_c(conditional|L); ll_i(integrated); ll_ei (expected i)
        // --------------------------------------------------------------------------
        // ll.fill(0.0);
        // for(uword ak=0; ak<K; ak++){
        //   ll.row(ak) = epmodel.log_p(epdta, ak);
        // }
        // ll_c += accu(par.L(sub).slice(ie) % ll2) / dta.np;
        // ll_i += epmodel.avg_log_p(epdta);
        // epmodel.set_params(post.Emeans, post.Edcovs, mpij);
        // ll_ei+= epmodel.avg_log_p(epdta);
      }
      // ----------------------------------------------------------------------------
      // SUBmodule 1: Ci | alpha, pi
      // ----------------------------------------------------------------------------
      mat llc = zeros(K, dta.np);
      // if(iter > 0) llc += nu_est(par.S, par.alpha(sub));         //MIC1
      if(iter > 0) llc.each_col() += arma::log(mpi.t());         //MIC2
      for(int ne=0; ne<dta.ne(sub); ne++) llc += nu_est(par.L(sub).slice(ne), par.beta(sub)(ne));
      par.C.slice(sub) = mrmultinom0(llc);    //use mrmultinom for MIC1
      if (iter > 0) clustalign(par.C.slice(sub), par.S);
      // ----------------------------------------------------------------------------
      // SUBmodule 2: Beta | C, L, pr
      // ----------------------------------------------------------------------------
      for(int ie2=0; ie2<dta.ne(sub); ie2++){
        par.beta(sub)(ie2) = sampcoh(par.L(sub).slice(ie2), par.C.slice(sub),
                 pr.a1, pr.b1);
      }
    }
    // ------------------------------------------------------------------------------
    // POPmodule 1: S | alpha, pi
    // ------------------------------------------------------------------------------
    mat lls = zeros(K, dta.np);               // ll mat for S
    for(int sp=0; sp<dta.np; sp++) lls.col(sp) += arma::log(par.pi);
    for(int ss=0; ss<dta.ns; ss++) lls += nu_est(par.C.slice(ss), par.alpha(ss));
    par.S = mrmultinom0(lls);
    // ----------------------------------------------------------------------------
    // POPmodule 2: Alpha | S, C, pr
    // ----------------------------------------------------------------------------
    for(int ie3=0; ie3<dta.ns; ie3++){
      par.alpha(ie3) = sampcoh(par.C.slice(ie3), par.S, pr.a1, pr.b1);
    }
    // ----------------------------------------------------------------------------
    // POPmodule 7: pi |
    // ----------------------------------------------------------------------------
    par.pi = rDir(sum(par.S, 1), pr.dir0[0]);
    //
    // // IC's
    // par.ll_c = ll_c / accu(dta.ne);
    // par.ll_i = ll_i / accu(dta.ne);
    // par.ll_ei= ll_ei/ accu(dta.ne);
    //
    //-----------------------------------------------------------------------------
    // Output samples with 1/5 burn-in
    //-----------------------------------------------------------------------------
    //
    // 1/5 burning in + thining
    if((iter>=run/5.0) && (iter%thin==0)){
      OutputSample();
      //
      // Epoch output
      for(int j=0; j<dta.ns; j++){
        std::stringstream ss;
        for(int ne=0; ne<dta.ne(j); ne++){
          for(int np=0; np<dta.np; np++){
            int outk=0; while(par.L(j).slice(ne)(outk, np) < 0.9) outk++;
            ss << (int) outk+1 << " ";
          }
          ss << par.beta(j)(ne) << "\n";
        }
        out3[j] << ss.str();
      }
      // output counter update;
      oc++;
    }
  }
  // Close output files -----------------------------------------------------------
  out1.close(); out2.close();
  for(int s = 0; s<dta.ns; s++) out3[s].close();
  delete [] out3;
  //
  // Summarize posteior samples ---------------------------------------------------
  //
  PostSummary(oc);                      // All par's updated as the posterior mean
  // Estimate ll_e based on the posterior L;
  {
    double ll_e = 0;
    for(int sub=0; sub<dta.ns; sub++){
      rowvec mpi(K);
      par.pi /= arma::accu(par.pi);
      mpi = mixpr(par.pi, par.alpha(sub));
      // --------------------------Epoch level---------------------------------------
      for(int ie=0; ie<dta.ne(sub); ie++){
        subprior spr;
        spr.mu0     = pr.mu0(sub).col(ie);      spr.dcov0 = pr.dcov0(sub).col(ie);
        spr.b0      = pr.b0;
        mat epdta = dta.dta(sub).slice(ie);
        vec tempi = arma::conv_to< vec >::from(mpi);
        rowvec mpij = mixpr(tempi, par.beta(sub)(ie));
        mpij /= arma::accu(mpij);             // normal mix probs
        NIG post = NIGpost(epdta, spr, par.L(sub).slice(ie));
        arma::gmm_diag epmodel;
        mpij[0] = 1-accu(mpij(span(1,K-1)));  // numerical check...
        epmodel.set_params(post.Emeans, post.Edcovs, mpij);
        mat ll(K, dta.np);
        for(int k=0; k<K; k++){
          ll.row(k) = epmodel.log_p(epdta, k);
        }
        ll_e += arma::accu(par.L(sub).slice(ie) % ll) / (dta.np+0.0);
      }
    }
    par.ll_e = ll_e / arma::accu(dta.ne);
  }
  // IC's calculation
  vec ICs(2,fill::zeros);
  {
    double total_e = arma::accu(dta.ne);
    double BIC_c = (2.0 * K * arma::accu(dta.d % dta.ne) + K - 1.0) *
      (std::log(dta.np+0.0) + std::log(total_e)) / total_e / (dta.np+0.0);
    ICs[0] = 2.0 * par.ll_e - BIC_c;      //BIC
    // ICs[1] = 2.0 * par.ll_ei - BIC_c;     //BIC2
    // ICs[2] = 2.0 * par.ll_i - par.ll_ei;  //DIC4
    // ICs[3] = 2.0 * par.ll_c - par.ll_e;   //DIC7
    ICs[1] = (K * mean(par.alpha) - 1.0) / (K - 1.0);
  }
  //
  // Returning --------------------------------------------------------------------
  return List::create(Named("alpha")  = par.alpha,
                      Named("beta")   = par.beta,
                      Named("pi")     = par.pi,
                      Named("clustS") = par.S,
                      Named("clustC") = par.C,
                      Named("clustL") = par.L,
                      Named("ICs")    = ICs
    );
}
