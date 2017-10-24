#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec inferPPA(vec pval,
             const mat& ann,
             const vec& ann_w,
             double ann_w0=-4.6,
             double alpha=0.1,
             double minPval=1e-30);

vec inferPPAByLoci(const vec& zscore,
                   const mat& ann, 
                   const vec& ann_w, 
                   const vec& ann_w0, 
                   const vec& s2,
                   const umat& locusCursor);
