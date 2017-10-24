#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

double eval_lld(vec zscore, mat ldmat, double s2 );

double eval_config(const vec& zscore, 
                   urowvec config, 
                   const mat& ldmat, double s2);

double eval_config_hash(const vec& zscore, 
                        uvec config, 
                        std::unordered_map<std::vector<bool>, double>& configHash,
                        const mat& ldmat, 
                        const mat& ann, const vec& ann_w, 
                        double ann_w0, double s2);

IntegerVector sampleConfigs_nbh(const vec& zscore,
                                IntegerVector cur_config,
                                std::unordered_map<std::vector<bool>, double> &configHash,
                                const mat& ldmat,
                                const mat& ann, 
                                const vec& ann_w, double ann_w0,
                                double s2);

umat sampleConfigs(const vec& pip, int sampleSize);

