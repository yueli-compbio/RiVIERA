#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <unordered_map>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

#include "basicFuns.h"
#include "nbh.h"

// [[Rcpp::export]]
double eval_lld(vec zscore, mat ldmat, double s2 ) {
  
  mat causalmat = eye(zscore.n_elem, zscore.n_elem) * s2;
  
  double lld_causal = dmvnrm(zscore, 0, ldmat + ldmat * causalmat * ldmat);
  
  double lld_null = dmvnrm( zscore, 0, ldmat);
  
  double bf;
  
  if(is_finite(lld_causal) && is_finite(lld_null)) {
    
    bf = lld_causal - lld_null; 
    
  } else {
    
    Rcout << endl << "**eval_lld:" << endl;
    Rcout << "zscore: " << zscore << endl;
    Rcout << "ldmat: " << ldmat << endl;
    Rcout << "s2: " << s2 << endl;
    Rcout << "lld_causal: " << lld_causal << endl;
    Rcout << "lld_null: " << lld_null << endl;
    
    throw std::runtime_error("lld_causal/lld_null is not finite!");
    
    bf = -INFINITY;
  }
  
  // bayes factor
  return bf;
}


// [[Rcpp::export]]
double eval_pri(const mat& ann, const uvec& causalidx, 
                const vec& ann_w, double ann_w0, bool logd = true) {
  
  vec config_pri = 1 - inferPri(ann, ann_w, ann_w0);
  
  config_pri(causalidx) = 1 - config_pri(causalidx);
  
  double res;
  
  if(logd) {
    res = accu(trunc_log(config_pri));
  } else {
    res = prod(config_pri);
  }
  
  return res;
}


double eval_config(const vec& zscore, 
                   urowvec config, 
                   const mat& ldmat, double s2)
{
  uvec causalIdx = find(config > 0);
  
  return eval_lld(zscore( causalIdx ), ldmat( causalIdx, causalIdx ), s2);
}


double eval_config_hash(const vec& zscore, 
                        uvec config, 
                        std::unordered_map<std::vector<bool>, double>& configHash,
                        const mat& ldmat, 
                        const mat& ann, const vec& ann_w, 
                        double ann_w0, double s2) 
{
  int m = zscore.n_elem;
  
  double config_val = 0;
  
  std::unordered_map<std::vector<bool>, double>::const_iterator it;
  
  std::vector<bool> mykey(m, false);
  
  for(int j = 0; j < (int)config.n_elem; j++) {
    
    mykey[config(j)] = true;
  }
  
  it = configHash.find(mykey);
  
  // already in hash table
  if(it != configHash.end()) {
    
    config_val = it->second;
    
  } else {
    
    config_val =
      eval_lld( zscore( config ),
                ldmat( config, config ), s2 ) +
                  eval_pri( ann, config, ann_w, ann_w0 );
    
    
    if(is_finite(config_val)) {
      
      configHash[mykey] = config_val;
    }
  }
  
  return config_val;
}




// [[Rcpp::export]]
umat find_uniqueRows(umat x) {
  
  std::unordered_map<std::vector<bool>, int> uniqConfigs;
  
  int m = x.n_cols;
  
  for(int i=0; i<(int)x.n_rows; i++) {
    
    std::vector<bool> mykey(m, false);
    
    for(int j = 0; j < m; j++) {
      
      mykey[j] = x(i,j);
    }
    
    // excluding all-zero row
    if(any(x.row(i) > 0)) {
      uniqConfigs[mykey] = i; 
    }
  }
  
  std::unordered_map<std::vector<bool>, int>::iterator
    lhs = uniqConfigs.begin(), rhs = uniqConfigs.end();
  
  std::unordered_map<std::vector<bool>, int>::const_iterator it;
  
  uvec uniqueRowIdx(uniqConfigs.size());
  
  for (R_xlen_t i = 0; lhs != rhs; ++lhs, i++) {
    
    std::vector<bool> key_i = lhs->first;
    
    it = uniqConfigs.find(key_i);
    
    uniqueRowIdx(i) = it->second;
  }
  
  
  return x.rows(uniqueRowIdx);
}


// [[Rcpp::export]]
umat sampleConfigs(const vec& proposal, int sampleSize)
{
  // output dim: sampled configs x snps
  return find_uniqueRows(repmat(proposal.t(),sampleSize, 1) > 
                           randu(sampleSize, proposal.n_elem));
}



IntegerVector sampleConfigs_nbh(const vec& zscore,
                                IntegerVector cur_config,
                                std::unordered_map<std::vector<bool>, double> &configHash,
                                const mat& ldmat,
                                const mat& ann, 
                                const vec& ann_w, double ann_w0,
                                double s2)
{
  int m = zscore.n_elem;
  
  // generate neighborhood of the current config
  List cur_config_nbh = nbh(cur_config, m);
  
  int nbhSize = cur_config_nbh.size();
  
  NumericVector cur_nbh_prob(nbhSize);
  
  for(int k=0; k < nbhSize; k++) {
    
    cur_nbh_prob[k] = eval_config_hash(zscore, as<uvec>(cur_config_nbh[k]), 
                                       configHash,
                                       ldmat, ann, ann_w, ann_w0, s2);
  }
  
  // convert to zero-based
  IntegerVector configid = seq_len(nbhSize) - 1;
  
  
  // divide by the largest number of all prob
  // before sum to get normalization to avoid inf
  cur_nbh_prob = exp(cur_nbh_prob - max(cur_nbh_prob));
  
  cur_nbh_prob = cur_nbh_prob / sum(cur_nbh_prob);
  
  if(any(is_na(cur_nbh_prob))) {
    for (int k = 0; k < nbhSize; k++) {
      if(NumericVector::is_na(cur_nbh_prob[k])) {
        cur_nbh_prob[k] = 0;  
      }
    } 
  }
  
  int new_idx = 0;
  
  try{
    
    new_idx = as<int>(Rcpp::RcppArmadillo::sample(
      
      configid, 1, false, cur_nbh_prob));
    
  } catch(...) {
    
    new_idx = as<int>(Rcpp::RcppArmadillo::sample(
      
      configid, 1, false));
  }
  
  return cur_config_nbh[new_idx];
}


































