#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec hmcS(const vec& zscore, 
         const List& ldmat_list,
         const List& configs_list, 
         const List& configs_post_list,
         const umat& locusCursor,
         vec s2_cur,
         double step, int nsteps,
         bool verbose, 
         double lbs2 = 1e-8, 
         double ubs2 = 1e3);
