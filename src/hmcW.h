#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec hmcW(const mat& ann,
         const List& configs_list, 
         const List& configs_post_list,
         const umat& locusCursor,
         vec ann_w_cur,
         vec ann_w_mu,
         mat ann_w_lambda,
         vec ann_w0_cur,
         vec ann_w0_mu,
         double ann_w0_tau,
         double step, int nsteps, 
         bool verbose);
