#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List inferConfigs_imp(const vec& zscore, 
                      const mat& ldmat, 
                      const mat& ann, 
                      const vec& ann_w, 
                      double ann_w0, 
                      double s2,
                      int sampleSize,
                      bool verbose);

List inferConfigs_nbh(const vec& zscore, 
                      const mat& ldmat, 
                      const mat& ann, 
                      const vec& ann_w, 
                      double ann_w0, 
                      double s2,
                      int sampleSize, 
                      bool verbose);

List inferConfigsAllLoci(const vec& proposal,
                         const vec& zscore, 
                         const List& ldmat_list, 
                         const umat& locusCursor,
                         const mat& ann, 
                         const vec& ann_w, 
                         const vec& ann_w0, 
                         const vec& s2,
                         int sampleSize,
                         int sampleScheme,
                         bool verbose=false);

vec inferPIP(const List& configs_list,
             const List& configs_post_list,
             const umat& locusCursor,
             int m, int locNum);

