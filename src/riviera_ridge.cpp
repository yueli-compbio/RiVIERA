#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "basicFuns.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec calGradPT_ridge(const vec& pri,
                    const vec& ppa,
                    const mat& ann,
                    const vec& ann_w,
                    const vec& ann_w_mu,
                    double ann_w0,
                    double ann_w0_mu)
{
  // vec pri = inferPri(ann, ann_w, ann_w0);
  
  int K = ann.n_cols;
  
  vec grad(K+1);
  
  // intercept
  grad(0) = accu(ppa - pri) - (ann_w0 - ann_w0_mu);
  
  // weights
  grad(span(1,K)) = ann.t() * (ppa - pri) - (ann_w - ann_w_mu);
  
  return grad;
}

// [[Rcpp::export]]
mat calHessPT_ridge(const vec& pri, const mat& ann)
{
  // vec pri = inferPri(ann, ann_w, ann_w0);
  
  int K = ann.n_cols;
  
  mat hess = zeros<mat>(K+1,K+1);
  
  hess(0,0) = as_scalar(pri.t() * (pri - 1));
  
  hess(span(1,K), span(0,0)) = ann.t() * (pri % (pri - 1));
  
  hess(span(0,0), span(1,K)) = hess(span(1,K), span(0,0)).t();
  
  hess(span(1,K),span(1,K)) = ann.t() * (ann % repmat(pri % (pri - 1), 1, K));
  
  hess.diag() = hess.diag() - 1;
  
  return hess;
}


// [[Rcpp::export]]
double lprPT_ridge(vec pval,
                   vec pri,
                   mat ann,
                   double alpha,
                   vec ann_w,
                   vec ann_w_mu,
                   double ann_w0,
                   double ann_w0_mu)
{
  // vec pri = inferPri(ann, ann_w, ann_w0);
  
  return accu((trunc_log(pri) + lbetaPT(pval, alpha)) + trunc_log(1-pri)) -
    accu(square(ann_w - ann_w_mu))/2 - pow(ann_w0 - ann_w0_mu,2)/2;
}


// [[Rcpp::export]]
double lprPT_givenPPA_ridge(vec pval,
                            vec pri,
                            vec ppa,
                            mat ann,
                            double alpha,
                            vec ann_w,
                            vec ann_w_mu,
                            double ann_w0,
                            double ann_w0_mu)
{
  // vec pri = inferPri(ann, ann_w, ann_w0);
  
  return accu( ppa % (trunc_log(pri) + lbetaPT(pval, alpha)) + (1-ppa) % trunc_log(1-pri) ) -
    accu(square(ann_w - ann_w_mu))/2 - pow(ann_w0 - ann_w0_mu,2)/2;
  
  // return
  //   accu(lbetaPT(pval, alpha)) +
  //     accu(ppa % trunc_log(pri) + (1-ppa) % trunc_log(1-pri)) -
  //     accu(square(ann_w - ann_w_mu))/2 - pow(ann_w0 - ann_w0_mu,2)/2;
}


// [[Rcpp::export]]
List riviera_ridge(vec pval,
                   const mat& ann,
                   vec ann_w_mu,
                   double minPval=1e-30,
                   double alpha0=0.1,
                   double pri0 = 0.1,
                   double epsilon=1e-5,
                   double thres=1e-3,
                   int max_iter=1e3,
                   double minAlpha=1e-30,
                   double maxAlpha=0.5,
                   bool fitAlpha=true)
{
  double alpha = alpha0;
  
  List emfit;
  
  vec lpr_init;
  
  uvec smallPval_idx = find(pval < minPval);
  
  int smallPval_cnt = smallPval_idx.n_elem;
  
  if(smallPval_cnt > 0) {
    
    Rcout << smallPval_cnt << 
      " SNPs have p-value < " << minPval << endl;
    
    Rcout << "They are set to " << minPval;
    Rcout << " for numerical stability" << endl;
    
    pval(smallPval_idx).fill(minPval);  
  }
  
  double ann_w0_mu = trunc_log(pri0) - trunc_log(1-pri0);
  
  double ann_w0 = ann_w0_mu;
  
  vec ann_w = ann_w_mu + 1e-5 * randn<vec>(ann.n_cols);
  
  // Full EM update
  Rcout << "Estimate parameters by model by EM:" << endl;
  
  // intialize
  int N = pval.n_elem;
  
  int iter = 0;
  
  double lpr_diff = 1;
  
  vec lpr = zeros<vec>(max_iter);
  
  // E-step
  vec pri = inferPri(ann, ann_w, ann_w0);
  
  vec causal = pri % dbetaPT(pval, alpha);
  
  vec ppa = causal/(causal + (1-pri));
  
  double N1 = accu(ppa); // causal SNPs count
  
  // Evaluate
  lpr(iter) = lprPT_ridge(pval,pri,ann,alpha,ann_w,ann_w_mu,ann_w0,ann_w0_mu);
  
  // lpr(iter) = lprPT_givenPPA_ridge(pval,pri,ppa,ann,alpha,ann_w,ann_w_mu,ann_w0,ann_w0_mu);
  
  Rcout << "iter " << iter+1 << ": lpr = " << lpr(iter);
  Rcout << "; N1 = " << N1;
  Rcout << "; alpha = " << alpha << "; pri0 = " << 1/(1+exp(-ann_w0));
  Rcout << "; diff = " << lpr_diff << endl;
  
  iter++;
  
  vec lpval = trunc_log(pval);
  
  int K = ann.n_cols;
  
  vec newtonUpdates(K+1);
  
  while(iter < max_iter && lpr_diff > thres) {
    
    // M-step
    if(fitAlpha) {
      
      alpha = -N1/accu(ppa % lpval);
      
      if(alpha > maxAlpha) {
        alpha = maxAlpha;
      } else if (alpha < minAlpha) {
        alpha = minAlpha;
      }
    }
    
    // gradient update
    // ann_w += epsilon * (ann.t() * (pri - ppa) + (ann_w - ann_w_mu));
    // ann_w0 += epsilon * (accu(pri - ppa) + (ann_w0 - ann_w0_mu));
    
    // newton update: inv(hessian) * gradient
    // vec newtonUpdates = pinv(hess) * grad;
    newtonUpdates = pinv(calHessPT_ridge(pri, ann)) * 
      calGradPT_ridge(pri, ppa, ann, ann_w, ann_w_mu, ann_w0, ann_w0_mu);
    
    ann_w0 -= epsilon * newtonUpdates(0);
    
    ann_w -= epsilon * newtonUpdates(span(1,K));
    
    pri = inferPri(ann, ann_w, ann_w0);
    
    // Evaluate
    lpr(iter) = lprPT_ridge(pval,pri,ann,alpha,ann_w,ann_w_mu,ann_w0,ann_w0_mu);
    
    // lpr(iter) = lprPT_givenPPA(pval,pri,ppa,ann,alpha,ann_w,ann_w_mu,ann_w0,ann_w0_mu);
    
    lpr_diff = (lpr(iter) - lpr(iter-1))/lpr(iter);
    
    if(lpr_diff<0) {
      lpr_diff = -lpr_diff;
    }
    
    Rcout << "iter " << iter+1 << ": lpr = " << lpr(iter);
    Rcout << "; N1 = " << N1;
    Rcout << "; alpha = " << alpha << "; pri0 = " << 1/(1+exp(-ann_w0));
    Rcout << "; diff = " << lpr_diff << endl;
    
    // E-step
    causal = pri % dbetaPT(pval, alpha);
    
    ppa = causal/(causal + (1-pri));
    
    N1 = accu(ppa); // causal SNPs count
    
    iter++;
  }
  
  lpr = lpr(span(0,iter-1));
  
  // covariance is the inverse of negative hessian matrix
  mat ann_w_cov = pinv(-calHessPT_ridge(inferPri(ann, ann_w, ann_w0), ann));
  
  // standard error
  vec ann_w_se = sqrt(abs(ann_w_cov.diag()));
  
  vec ann_w_zscore = ann_w/ann_w_se(span(1,K));
  
  return List::create(Named("alpha")=alpha,
                      Named("pri0")=pri0,
                      Named("ann_w")=ann_w,
                      Named("ann_w0")=ann_w0,
                      Named("ann_w_se")=ann_w_se,
                      Named("ann_w_cov")=ann_w_cov,
                      Named("ann_w_zscore")=ann_w_zscore,
                      Named("lpr_init")=lpr_init,
                      Named("lpr")=lpr);
}



















































