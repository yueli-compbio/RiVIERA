#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec calGradPT_ridge(const vec& pri,
                    const vec& ppa,
                    const mat& ann,
                    const vec& ann_w,
                    const vec& ann_w_mu,
                    double ann_w0,
                    double ann_w0_mu);

mat calHessPT_ridge(const vec& pri, const mat& ann);
