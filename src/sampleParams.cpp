#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "basicFuns.h"

// [[Rcpp::export]]
mat sampleAnnLambda(double wishart_df, mat ann_w_sigma0, vec ann_w) {
  
  return rwishart(wishart_df, (ann_w_sigma0 + ann_w * ann_w.t()).i());
}


// Gibbs sampling precision from gamma conjugate to
// the Gaussian distribution that bias follows
// double alpha=1e-2, double beta=1e-4
// [[Rcpp::export]]
double sampleAnnBiasTau(double ann_w0, double ann_w0_mu,
                        double alpha=1, double beta=1)
{
  return as<double>(rgamma(
      1, alpha + 1/2,
      1/(beta + pow( ann_w0 - ann_w0_mu, 2)/2)));
}
