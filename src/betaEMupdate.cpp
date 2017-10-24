#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "basicFuns.h"

using namespace Rcpp;
using namespace arma;


double lpr_estimateParByEM(vec pval, vec gamma, double alpha, double pri0) {
  
  // evaluate complete log likelihood
  return accu(trunc_log(pri0) + lbetaPT(pval, alpha)) + trunc_log(1-pri0);
}

// [[Rcpp::export]]
List estimateParByEM(vec pval,
                     double alpha_init = 0.1,
                     double pri0_init = 1e-3,
                     double minPval=1e-30,
                     double thres=1e-3,
                     int max_iter=1e3) 
{
  
  uvec smallPval_idx = find(pval < minPval);
  
  int smallPval_cnt = smallPval_idx.n_elem;
  
  if(smallPval_cnt > 0) {
    
    Rcout << smallPval_cnt << 
      " SNPs have p-value < " << minPval << endl;
    
    Rcout << "They are set to " << minPval;
    Rcout << " for numerical stability" << endl;
    
    pval(smallPval_idx).fill(minPval);  
  }
  
  double alpha = alpha_init;
  
  double pri0 = pri0_init;
  
  double lpr_diff = 1;
  
  int N = pval.n_elem;
  
  vec lpval = trunc_log(pval);
  
  int iter=0;
  
  vec lpr = zeros<vec>(max_iter);
  
  // E-step
  vec causal = pri0 * dbetaPT(pval, alpha);
  
  vec gamma = causal/(causal + (1-pri0));
  
  double N1 = accu(gamma); // causal SNPs count
  
  lpr(iter) = lpr_estimateParByEM(pval, gamma, alpha, pri0);
  
  Rcout << "iter " << iter+1 << ": lpr = " << lpr(iter);
  Rcout << "; alpha = " << alpha << "; pri0 = " << pri0;
  Rcout << "; diff = " << lpr_diff << endl;
  iter++;
  
  double alphaMIN = 1e-5;
  double alphaMAX = 1;
  
  while(iter < max_iter && lpr_diff > thres) {
    
    // M-step
    alpha = -N1/accu(gamma % lpval);
    
    if(alpha < alphaMIN) alpha = alphaMIN;
    if(alpha > alphaMAX) alpha = alphaMAX;
    
    pri0 = N1/N;
    
    lpr(iter) = lpr_estimateParByEM(pval, gamma, alpha, pri0);
    
    lpr_diff = lpr(iter) - lpr(iter-1);
    
    Rcout << "iter " << iter+1 << ": lpr = " << lpr[iter];
    Rcout << "; alpha = " << alpha << "; pri0 = " << pri0;
    Rcout << "; diff = " << lpr_diff << endl;
    
    // E-step
    causal = pri0 * dbetaPT(pval, alpha);
    
    gamma = causal/(causal + (1-pri0));
    
    N1 = accu(gamma); // causal SNPs count
    
    iter++;
  }
  
  lpr = lpr(span(0,iter-1));
  
  return List::create(Named("alpha")=alpha,
                      Named("pri0")=pri0,
                      Named("lpr")=lpr);
}