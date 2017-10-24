#include <RcppArmadillo.h>
#include <unordered_map>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "basicFuns.h"
#include "sampleConfigs.h"
#include "calGrad.h"


// [[Rcpp::export]]
List inferConfigs_imp(const vec& proposal,
                      const vec& zscore, 
                      const mat& ldmat, 
                      const mat& ann, 
                      const vec& ann_w, 
                      double ann_w0, 
                      double s2,
                      int sampleSize,
                      bool verbose)
{
  int m = zscore.n_elem;
  
  // evaluate ppa to form the sampling basis
  vec pri = inferPri(ann, ann_w, ann_w0);
  
  // bayesfactor
  // vec bfr = (pri % duvnrm(zscore, 0, 1+s2))/((1-pri) % duvnrm(zscore));
  // vec ppa = bfr/accu(bfr);
  
  // n x m
  umat configs = sampleConfigs(proposal, sampleSize);
  
  int n = configs.n_rows;
  
  if(verbose) {
    
    Rcout << n << " configs were sampled from ";
    Rcout  << m << " SNPs" << endl << endl;
  }
  
  vec config_lld(n);
  
  for(int i=0; i<n; i++) {
    
    config_lld(i) = eval_config(zscore, configs.row(i), ldmat, s2);
  }
  
  // offset sampling bias towards snps with large bayesfactor
  // log likelihood + log importance weight
  vec configs_post = config_lld +
    configs * (trunc_log(pri) - trunc_log(proposal)) +
    (1-configs) * (trunc_log(1-pri) - trunc_log(1-proposal));
  
  // vec configs_post = config_lld + 
  //   configs * trunc_log(pri) + (1-configs) * trunc_log(1-pri);
  
  // divide by max to avoid inf sum
  double configVal_max = max(configs_post);
  
  configs_post = configs_post - configVal_max;
  
  configs_post = trunc_exp(configs_post);
  
  configs_post = configs_post / sum(configs_post);
  
  List res = List::create(
    Named("configs") = configs,
    Named("configs_post") = configs_post);
  
  return res;
}


// [[Rcpp::export]]
List inferConfigs_nbh(const vec& zscore, 
                       const mat& ldmat, 
                       const mat& ann, 
                       const vec& ann_w, 
                       double ann_w0, 
                       double s2,
                       int sampleSize, 
                       bool verbose)
{
  int m = zscore.n_elem;
  
  // hash table to record all sampled valid configs
  std::unordered_map<std::vector<bool>, double> configHash;
  
  IntegerVector config(0);
  
  int configHashSize = 0;
  
  for(int t=0; t < sampleSize; t++) {
    
    config = sampleConfigs_nbh(zscore,
                               config,
                               configHash,
                               ldmat, ann,
                               ann_w, ann_w0, s2);
    
    configHashSize = configHash.size();
  }
  
  urowvec myconfig = zeros<urowvec>(m);
  
  configHashSize = configHash.size();
  
  if(verbose) {
    Rcout << configHashSize << " configs were sampled from ";
    Rcout << m << " SNPs" << endl;
  }
  
  std::unordered_map<std::vector<bool>, double>::iterator
    lhs = configHash.begin(), rhs = configHash.end();
  
  // save for later retrieval when sampling s2
  umat configs(configHashSize, m);
  
  // save for later retrieval when sampling s2
  vec configs_post(configHashSize);
  
  for (R_xlen_t i = 0; lhs != rhs; ++lhs, i++) {
    
    std::vector<bool> key_i = lhs->first;
    
    for(int j=0; j<m; j++) {
      
      myconfig(j) = key_i[j];
    }
    
    configs.row(i) = myconfig;
    
    configs_post(i) = lhs->second;
  }
  
  configHash.clear();
  
  // divide by max to avoid inf sum
  configs_post = trunc_exp(configs_post - max(configs_post));
  
  configs_post = configs_post / accu(configs_post);
  
  return List::create(
    Named("configs") = configs,
    Named("configs_post") = configs_post);
}



// [[Rcpp::export]]
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
                         bool verbose=false) 
{
  // target global cursor across all snps over all loci
  int global_up = 0;
  int global_dw = 0;
  
  int locNum = ldmat_list.size();
  
  double lprS = 0;
  
  List res;
  List configs_list;
  List configs_post_list;
  
  int totalConfigs = 0;
  
  for(int l=0; l < locNum; l++) {
    
    global_up = locusCursor(l, 0);
    
    global_dw = locusCursor(l, 1);

    if(sampleScheme==1) {
      
      res = inferConfigs_nbh(zscore(span(global_up, global_dw)), 
                             as<mat>(ldmat_list[l]), 
                             ann.rows(span(global_up, global_dw)), 
                             ann_w, 
                             ann_w0(l), 
                             s2(l),
                             sampleSize,
                             verbose);
    } else {
      
      res = inferConfigs_imp(proposal(span(global_up, global_dw)),
                             zscore(span(global_up, global_dw)), 
                             as<mat>(ldmat_list[l]), 
                             ann.rows(span(global_up, global_dw)), 
                             ann_w, 
                             ann_w0(l), 
                             s2(l),
                             sampleSize,
                             verbose);
    }
    
    configs_list.push_back(as<umat>(res["configs"]));
    configs_post_list.push_back(as<vec>(res["configs_post"]));
  }
  
  return List::create(
    Named("configs_list") = configs_list,
    Named("configs_post_list") = configs_post_list );
}



vec inferPIP(const List& configs_list,
             const List& configs_post_list,
             const umat& locusCursor,
             int m, int locNum) 
{
  int global_up = 0;
  int global_dw = 0;
  
  vec pip(m);
  
  for(int l=0; l<locNum; l++) {
    
    global_up = locusCursor(l, 0);
    
    global_dw = locusCursor(l, 1);
    
    pip(span(global_up, global_dw)) = 
      
      as<umat>(configs_list[l]).t() * as<vec>(configs_post_list[l]);
  }
  
  return pip;
}















