#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "basicFuns.h"
#include "sampleConfigs.h"

// log prior of model params
double priW(const vec& ann_w, 
            const vec& ann_w_mu,
            const mat& ann_w_lambda,
            
            vec ann_w0,
            vec ann_w0_mu,
            double ann_w0_tau) 
{
  
  // ann_w prior
  double prior_ann_w = -as_scalar((ann_w.t() - ann_w_mu.t()) *
                                  ann_w_lambda * (ann_w - ann_w_mu))/2;
  
  double prior_ann_w0 = accu(-ann_w0_tau * pow(ann_w0 - ann_w0_mu, 2)/2);
  
  
  // log gaussian prior + log causal likelihood
  return prior_ann_w + prior_ann_w0;
}

// log prior of model params
double priW_l(const vec& ann_w, 
              const vec& ann_w_mu,
              const mat& ann_w_lambda,
              
              double ann_w0,
              double ann_w0_mu,
              double ann_w0_tau) 
{
  
  // ann_w prior
  double prior_ann_w = -as_scalar((ann_w.t() - ann_w_mu.t()) *
                                  ann_w_lambda * (ann_w - ann_w_mu))/2;
  
  double prior_ann_w0 = -ann_w0_tau * pow(ann_w0 - ann_w0_mu, 2)/2;
  
  // log gaussian prior + log causal likelihood
  return prior_ann_w + prior_ann_w0;
}


// log posterior of model params
double likW_l(const vec& pri,
              const umat& configs,
              const vec& configs_post)
{
  return as_scalar(
    configs_post.t() * 
      (configs * trunc_log(pri) + (1-configs) * trunc_log(1-pri))
  );
}

double lprW_l(const vec& pri,
              const umat& configs,
              const vec& configs_post,
              
              const vec& ann_w,
              const vec& ann_w_mu,
              const mat& ann_w_lambda,
              
              double ann_w0,
              double ann_w0_mu,
              double ann_w0_tau)
{
  return likW_l(pri, configs, configs_post) + 
    priW_l(ann_w, ann_w_mu, ann_w_lambda, ann_w0, ann_w0_mu, ann_w0_tau);
}


// log posterior of model params
// [[Rcpp::export]]
double lprW(const mat& ann,
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
  double lld = 0;
  
  int global_up = 0;
  int global_dw = 0;
  
  for(int l=0; l<(int)configs_post_list.size(); l++) {
    
    global_up = locusCursor(l, 0);
    global_dw = locusCursor(l, 1);
    
    // log causal likelihood
    lld += likW_l(inferPri(ann.rows(span(global_up, global_dw)), 
                           ann_w, ann_w0(l)), 
                  as<umat>(configs_list[l]), 
                  as<vec>(configs_post_list[l]));
  }
  
  return lld + priW(ann_w, ann_w_mu, ann_w_lambda, 
                    ann_w0, ann_w0_mu, ann_w0_tau);
}


// s2 terms of the lpr
// [[Rcpp::export]]
double lprS_l(const vec& zscore,
              const mat& ldmat,
              const umat& configs,
              const vec configs_post,
              double s2)
{
  double lprS = 0;
  
  int configNum = configs_post.n_elem;
  
  for(int j = 0; j < configNum; j++) {
    
    uvec causalidx = find(configs.row(j)==1);
    
    // log expected complete posterior of causal configs
    lprS += configs_post(j) *
      eval_lld( zscore( causalidx ), ldmat( causalidx, causalidx ), s2 );
    
  }
  
  return lprS;
}



double lprS(const vec& zscore,
            const List& ldmat_list,
            const List& configs_list,
            const List& configs_post_list,
            const umat& locusCursor,
            const vec& s2)
{
  double lprS = 0;
  
  int global_up = 0;
  int global_dw = 0;
  
  for(int l=0; l<(int)configs_post_list.size(); l++) {
    
    global_up = locusCursor(l, 0);
    global_dw = locusCursor(l, 1);
    
    // log causal likelihood
    lprS += lprS_l(zscore(span(global_up, global_dw)),
                   as<mat>(ldmat_list[l]),
                   as<umat>(configs_list[l]),
                   as<vec>(configs_post_list[l]),
                   s2(l));
  }
  
  return lprS;
}






















