#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec calGradW(const mat& ann, 
             
             const List& configs_list, 
             const List& configs_post_list,
             const umat& locusCursor,
             
             const vec& ann_w,
             const vec& ann_w_mu,
             const mat& ann_w_lambda,
             
             vec ann_w0,
             vec ann_w0_mu,
             double ann_w0_tau);


vec calGradS(const vec& zscore, 
             const List& ldmat_list,
             const List& configs_list, 
             const List& configs_post_list,
             const umat& locusCursor,
             const vec& s2);
