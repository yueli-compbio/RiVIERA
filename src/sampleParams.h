#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

mat sampleAnnLambda(double wishart_df, mat ann_w_sigma0, vec ann_w);

double sampleAnnBiasTau(double ann_w0, double ann_w0_mu,
                        double alpha=1, double beta=1);
