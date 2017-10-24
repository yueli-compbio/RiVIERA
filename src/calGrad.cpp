#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "basicFuns.h"

// internal backprop function used by both backprop and hmc
// [[Rcpp::export]]
vec calGradW_l(const mat& ann, 
               
               const umat& configs, 
               const vec& configs_post,
               
               const vec& ann_w,
               const vec& ann_w_mu,
               const mat& ann_w_lambda,
               
               double ann_w0,
               double ann_w0_mu,
               double ann_w0_tau)
{
  int m = ann.n_rows;
  
  int K = ann.n_cols;
  
  // intercept + annotations
  vec grad = zeros<vec>(K + 1);
  
  // gradients for logistic coef
  // prior
  grad(0) = (ann_w0 - ann_w0_mu) * ann_w0_tau;
  
  grad(span(1,K)) = ann_w_lambda * (ann_w - ann_w_mu);
  
  int configNum = configs_post.n_elem;
  
  vec pri = inferPri(ann, ann_w, ann_w0);
  
  // mxn nx1 -> mx1
  mat tmp = (repmat(pri, 1, configNum) - configs.t()) * configs_post;
  
  grad(0) += accu(tmp);
  
  // kxm mx1 -> kx1
  grad(span(1,K)) += ann.t() * tmp;
  
  return grad;
}

// [[Rcpp::export]]
vec calGradW(const mat& ann, 
             
             const List& configs_list, 
             const List& configs_post_list,
             const umat& locusCursor,
             
             const vec& ann_w,
             const vec& ann_w_mu,
             const mat& ann_w_lambda,
             
             vec ann_w0,
             vec ann_w0_mu,
             double ann_w0_tau)
{
  int m = ann.n_rows;
  
  int K = ann.n_cols;
  
  int locNum = locusCursor.n_rows;
  
  // intercepts per locus + annotations
  vec grad = zeros<vec>(locNum + K);
  
  // gradients for logistic coef
  // prior
  grad(span(0,locNum-1)) = (ann_w0 - ann_w0_mu) * ann_w0_tau;
  
  grad(span(locNum, locNum+K-1)) = ann_w_lambda * (ann_w - ann_w_mu);
  
  // target global cursor across all snps over all loci
  int global_up = 0;
  int global_dw = 0;
  
  // Rcout << "pass" << endl;
  
  for(int l=0; l < locNum; l++) {
    
    // Rcout << "loc: " << l << endl;
    
    global_up = locusCursor(l, 0);
    global_dw = locusCursor(l, 1);
    
    int configNum = as<vec>(configs_post_list[l]).n_elem;
    
    vec pri = inferPri(ann.rows(span(global_up, global_dw)), ann_w, ann_w0(l));
    
    // mxn nx1 -> mx1
    mat tmp = (repmat(pri, 1, configNum) - as<umat>(configs_list[l]).t()) * 
      
      as<vec>(configs_post_list[l]);
    
    grad(l) += accu(tmp);
    
    // kxm mx1 -> kx1
    grad(span(locNum, locNum+K-1)) += ann.rows(span(global_up, global_dw)).t() * tmp;
  }
  
  return grad;
}



// [[Rcpp::export]]
mat calGrad_dEdV(vec x, mat sigma) {
  
  // mat rooti = trans(inv(trimatu(chol(sigma))));
  
  mat sigma_inv(size(sigma));
  
  try {
    
    sigma_inv = fastInv(sigma);
    
  } catch(...) {
    
    try {
      
      sigma_inv = pinv(sigma);
      
    } catch(...) {
      
      throw std::runtime_error("calGrad_dEdV: sigma is not invertible!");
    }
  }
  
  return 
    
    -0.5 * (sigma_inv - (sigma_inv * x * trans(x) * sigma_inv));
}


// [[Rcpp::export]]
double calGrad_dEdS(vec x, mat sigma, double s) {
  
  // dE/dS = dE/dV dV/dS
  return trace(calGrad_dEdV(x, sigma + 
               
               sigma * (eye(sigma.n_rows, sigma.n_rows) * s) * sigma) * 
               
               sigma * sigma);
}


// [[Rcpp::export]]
double calGradS_l(const vec& zscore, 
                  const mat& ldmat,
                  const umat& configs, 
                  const vec& configs_post,
                  double s2) 
{
  int configNum = configs_post.n_elem;
  
  double grad = 0;
  
  for(int j = 0; j < configNum; j++) {
    
    uvec causalidx = find(configs.row(j)==1);
    
    grad -= configs_post(j) * 
      calGrad_dEdS( zscore(causalidx), ldmat(causalidx, causalidx), s2 );
  }
  
  return grad;
}

// [[Rcpp::export]]
vec calGradS(const vec& zscore, 
             const List& ldmat_list,
             const List& configs_list, 
             const List& configs_post_list,
             const umat& locusCursor,
             const vec& s2) 
{
  double lprS = 0;
  
  int global_up = 0;
  int global_dw = 0;
  
  vec grad = zeros<vec>(s2.n_elem);
  
  for(int l=0; l<(int)locusCursor.n_rows; l++) {
    
    global_up = locusCursor(l, 0);
    global_dw = locusCursor(l, 1);
    
    // log causal likelihood
    grad(l) += calGradS_l(
      zscore(span(global_up, global_dw)),
      as<mat>(ldmat_list[l]),
      as<umat>(configs_list[l]),
      as<vec>(configs_post_list[l]),
      s2(l));
  }
  
  return grad;
}































