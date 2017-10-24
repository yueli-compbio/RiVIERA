#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

List nbh(IntegerVector config, int m);
