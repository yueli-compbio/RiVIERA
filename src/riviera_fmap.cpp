#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "basicFuns.h"
#include "hmcS.h"
#include "hmcW.h"
#include "inferConfigs.h"
#include "lpr.h"
#include "inferPPA.h"

// [[Rcpp::export]]
List riviera_fmap(
    const vec& zscore,
    const List& ldmat_list,
    const mat& ann, 
    const vec& ann_w_mu,
    const mat& ann_w_sigma0,
    int configSampleScheme=1,
    double c0=1.0,
    double pri0_init=1e-2,
    double s2_init=1,
    double burninFrac = 0.2,
    int sampleSize=1e3,
    int maxConfig = 1e6,
    int max_iter = 1e3,
    double step_W = 0.01, 
    int nsteps_W = 100,
    double step_S = 0.01, 
    int nsteps_S = 100,
    double thres = 1e-3,
    bool verbose = true,
    bool verboseC = true,
    bool verboseS = true,
    bool verboseW = true)
{
  Rcout << "Settings: " << endl;
  Rcout << "max_iter = " << max_iter << endl;
  Rcout << "step_W = " << step_W << "; nsteps_W = " << nsteps_W << endl;
  Rcout << "step_S = " << step_S << "; nsteps_S = " << nsteps_S << endl;
  Rcout << "configSampleScheme (1: nbh; 2: imp)= " << configSampleScheme << endl;
  
  int m = zscore.n_elem;
  int K = ann.n_cols;
  int locNum = ldmat_list.size();
  
  Rcout << m << " SNPs" << endl;
  Rcout << K << " annotations" << endl;
  Rcout << locNum << " loci" << endl;
  
  umat locusCursor = createLocusCursor(ldmat_list, m);
  
  vec pri0(locNum);

  for(int l=0; l<locNum; l++) {
    
    double m_l = locusCursor(l,1) - locusCursor(l,0) + 1.0;
    
    pri0(l) = c0/m_l;
  }

  pri0(find(pri0 > pri0_init)).fill(pri0_init);

  vec ann_w0_mu = log(pri0) - log(1-pri0);

  vec ann_w0 = ann_w0_mu;

  double ann_w0_tau = 1;

  // initialize enrichments to the prior mean
  vec ann_w = ann_w_mu;

  // precision matrix for w_k across K annotations
  mat ann_w_lambda = pinv(ann_w_sigma0);

  vec s2(locNum);
  s2.fill(s2_init);

  vec lpr_w = zeros<vec>(max_iter);
  vec accFreqW = zeros<vec>(max_iter);
  vec accRateW = zeros<vec>(max_iter);

  vec lpr_s = zeros<vec>(max_iter);
  vec accFreqS = zeros<vec>(max_iter);
  vec accRateS = zeros<vec>(max_iter);

  bool acceptFlagS = false; // indicate HMC acceptance
  bool acceptFlagW = false; // indicate HMC acceptance
  int ensembleSize = 0; // number of accepted models

  // save accepted models during the MCMC
  mat ann_w_ensemble = zeros<mat>(K, max_iter);

  mat ann_w0_ensemble = zeros<mat>(locNum, max_iter);

  mat s2_ensemble = ones<mat>(locNum, max_iter);

  mat pip_ensemble = zeros<mat>(m, max_iter);
  
  vec proposal = inferPPAByLoci(zscore, ann, ann_w, ann_w0, s2, locusCursor);
  
  // vec proposal = normalizeZscoreByLoci(zscore, locusCursor);
  
  // first sample use ppa instead of pip
  // E-step
  List res = inferConfigsAllLoci(proposal,
                                 zscore,
                                 ldmat_list,
                                 locusCursor,
                                 ann,
                                 ann_w,
                                 ann_w0,
                                 s2,
                                 sampleSize,
                                 configSampleScheme,
                                 verboseC);

  List configs_list = res["configs_list"];

  List configs_post_list = res["configs_post_list"];

  ann_w_ensemble.col(ensembleSize) = ann_w;

  ann_w0_ensemble.col(ensembleSize) = ann_w0;

  s2_ensemble.col(ensembleSize) = s2;

  pip_ensemble.col(ensembleSize) =
    inferPIP(configs_list, configs_post_list, locusCursor, m, locNum);

  ensembleSize++;

  vec ann_w_prev = ann_w;

  vec s2_prev = s2;

  vec retW(K+1);

  int iter = 0;
  
  // evaluate
  lpr_w(iter) = lprW(
    ann,
    configs_list,
    configs_post_list,
    locusCursor,
    ann_w, ann_w_mu, ann_w_lambda,
    ann_w0, ann_w0_mu, ann_w0_tau);
  

  lpr_s(iter) = lprS(zscore, ldmat_list,
        configs_list, configs_post_list, locusCursor, s2);
  

  double lprS_diff = 1;
  double lprW_diff = 1;
  double ann_w_diff = 1;

  // first sampling is special as model makes dramatic leaps
  if(verbose) {

    Rcout << iter+1;
    Rcout << ", lprS: " << lpr_s(iter);
    Rcout << "; lprS_diff: " << lprS_diff;
    Rcout << "; AFS: " << accFreqS(iter) << " (";
    Rcout << round(accRateS(iter)) << "%)";

    Rcout << ", lprW: " << lpr_w(iter);
    Rcout << "; lprW_diff: " << lprW_diff;

    Rcout << "; ann_w_diff: " << ann_w_diff;

    Rcout << "; AFW: " << accFreqW(iter) << " (";
    Rcout << round(accRateW(iter)) << "%)" << endl;
  }

  bool accepted = false;

  iter++;

  // keep track of snps at different loci
  int global_up = 0;
  int global_dw = 0;

  // infer marginal pip per locus per trait
  while(iter < max_iter && ann_w_diff > thres) {
    
    // Rcout << "Begin iter:" << iter+1 << endl;
    // Rcout << "configs_post_list: " << as<vec>(configs_post_list[0]) << endl;

    // M-step
    ann_w_prev = ann_w;
    
    s2_prev = s2;

    //// sample variance scaling parameters ////
    s2 = hmcS(zscore,
              ldmat_list,
              configs_list,
              configs_post_list,
              locusCursor,
              s2,
              step_S, nsteps_S, verboseS);
    
    // Rcout << "riviera_fmap: after hmcS, s2: " << s2(0) << endl;
    
    // HMC sample ann_w and ann_w0
    retW = hmcW(ann,
                configs_list,
                configs_post_list,
                locusCursor,
                ann_w,
                ann_w_mu,
                ann_w_lambda,
                ann_w0,
                ann_w0_mu,
                ann_w0_tau,
                step_W, nsteps_W, verboseW);
    
    ann_w0 = retW(span(0,locNum-1));
    
    ann_w = retW(span(locNum, locNum+K-1));

    acceptFlagS = any(s2_prev != s2);

    // collect fit info for diagnostic purpose
    if(iter > 0) {
      accFreqS(iter) = accFreqS(iter-1) + acceptFlagS;
    } else {
      accFreqS(iter) = acceptFlagS;
    }

    accRateS(iter) = 100*accFreqS(iter)/(iter+1);

    //// sample annotation weights ////
    // Gibb sample ann_w_lambda and ann_w0_tau
    // ann_w_lambda = sampleAnnLambda(ann_w_sigma0_estimateSize, ann_w_sigma0, ann_w);
    // ann_w0_tau = sampleAnnBiasTau(ann_w0, ann_w0_mu);

    acceptFlagW = any(ann_w!=ann_w_prev);

    // collect fit info for diagnostic purpose
    if(iter > 0) {
      accFreqW(iter) = accFreqW(iter-1) + acceptFlagW;
    } else {
      accFreqW(iter) = acceptFlagW;
    }

    accRateW(iter) = 100*accFreqW(iter)/(iter+1);

    // evaluate
    lpr_w(iter) = lprW(
      ann,
      configs_list,
      configs_post_list,
      locusCursor,
      ann_w, ann_w_mu, ann_w_lambda,
      ann_w0, ann_w0_mu, ann_w0_tau);

    lpr_s(iter) = lprS(zscore,
          ldmat_list,
          configs_list,
          configs_post_list,
          locusCursor, s2);

    lprW_diff = lpr_w(iter) - lpr_w(iter-1);

    lprS_diff = lpr_s(iter) - lpr_s(iter-1);

    accepted = acceptFlagS && acceptFlagW;

    // E-step
    if(accepted) {
      
      // Rcout << "riviera_fmap: If accepted, " << "s2: " << s2 << endl;

      s2_ensemble.col(ensembleSize) = s2;

      ann_w_ensemble.col(ensembleSize) = ann_w;

      ann_w0_ensemble.col(ensembleSize) = ann_w0;
      
      proposal = inferPPAByLoci(zscore, ann, ann_w, ann_w0, s2, locusCursor);
      
      res = inferConfigsAllLoci(proposal,
                                zscore,
                                ldmat_list,
                                locusCursor,
                                ann,
                                ann_w,
                                ann_w0,
                                s2,
                                sampleSize,
                                configSampleScheme,
                                verboseC);

      configs_list = res["configs_list"];

      configs_post_list = res["configs_post_list"];

      pip_ensemble.col(ensembleSize) =
        inferPIP(configs_list, configs_post_list, locusCursor, m, locNum);

      ann_w_diff = sqrt(accu(square(ann_w_ensemble.col(ensembleSize) -
        ann_w_ensemble.col(ensembleSize-1))))/
          sqrt(accu(square(ann_w_ensemble.col(ensembleSize) +
            ann_w_ensemble.col(ensembleSize-1))));

      ensembleSize++;
    }


    if(verbose) {

      Rcout << iter+1;
      Rcout << ", lprS: " << lpr_s(iter);
      Rcout << "; lprS_diff: " << lprS_diff;
      Rcout << "; AFS: " << accFreqS(iter) << " (";
      Rcout << round(accRateS(iter)) << "%)";

      Rcout << ", lprW: " << lpr_w(iter);
      Rcout << "; lprW_diff: " << lprW_diff;

      Rcout << "; ann_w_diff: " << ann_w_diff;

      Rcout << "; AFW: " << accFreqW(iter) << " (";
      Rcout << round(accRateW(iter)) << "%)" << endl;
    }

    iter++;
  }

  Rcout << "MCMC inference completed." << endl;

  Rcout << "Accepted models: " << ensembleSize << endl;

  pip_ensemble = pip_ensemble.cols(span(0,ensembleSize-1));
  ann_w_ensemble = ann_w_ensemble.cols(span(0,ensembleSize-1));
  ann_w0_ensemble = ann_w0_ensemble.cols(span(0,ensembleSize-1));
  s2_ensemble = s2_ensemble.cols(span(0,ensembleSize-1));

  lpr_w = lpr_w(span(0,iter-1));
  lpr_s = lpr_w(span(0,iter-1));

  List fit_info =
    List::create(
      Named("lprW")=lpr_w,
      Named("accFreqW")=accFreqW,
      Named("accRateW")=accRateW,
      Named("lprS")=lpr_s,
      Named("accFreqS")=accFreqS,
      Named("accRateS")=accRateS,
      Named("ensembleSize")=ensembleSize);


  List ensemble = List::create(
    Named("pip_ensemble")=pip_ensemble,
    Named("ann_w_ensemble")=ann_w_ensemble,
    Named("ann_w0_ensemble")=ann_w0_ensemble,
    Named("s2_ensemble")=s2_ensemble);

  return List::create(
    Named("ensemble")=ensemble,
    Named("fit_info")=fit_info);
}







