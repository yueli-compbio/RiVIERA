#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "basicFuns.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec inferPPA(vec pval,
             const mat& ann,
             const vec& ann_w,
             double ann_w0=-4.6,
             double alpha=0.1,
             double minPval=1e-30)
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
  
  vec causal = dbetaPT(pval, alpha);
  
  vec ppa(pval.n_elem);
  
  vec pri = inferPri(ann, ann_w, ann_w0);
  
  causal = pri % causal;
  
  ppa = causal/(causal + (1-pri));
  
  return ppa;
}

// [[Rcpp::export]]
vec inferPPAByLoci(const vec& zscore,
                   const mat& ann, 
                   const vec& ann_w, 
                   const vec& ann_w0, 
                   const vec& s2,
                   const umat& locusCursor)
{
  int global_up = 0;
  
  int global_dw = 0;
  
  vec ppa(ann.n_rows);
  
  for(int l=0; l<(int)locusCursor.n_rows; l++) {
    
    global_up = locusCursor(l, 0);
    
    global_dw = locusCursor(l, 1);
    
    vec pri = inferPri(ann.rows(span(global_up, global_dw)), ann_w, ann_w0(l));
    
    vec causal = pri % duvnrm(zscore(span(global_up, global_dw)), 0, 1+s2(l));
    
    vec null = (1-pri) % duvnrm(zscore(span(global_up, global_dw)));
    
    vec bfr = causal/null;
    
    ppa(span(global_up, global_dw)) = bfr/accu(bfr);
  }
  
  return ppa;
}





