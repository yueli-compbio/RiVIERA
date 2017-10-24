#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

List initParams(const List& zscoreList, 
                const List& annList, 
                const List& ldmatList, 
                const List& locusCursorList,
                double causalCntPerLocusPrior,
                double s2_init,
                double pri0=1e-2);
