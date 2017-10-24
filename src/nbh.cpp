#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


// add a causal causal to existing configuration
IntegerVector add(IntegerVector config, int snpidx, int causalcnt) {
  
  IntegerVector config_new(causalcnt+1);
  
  config_new[0] = snpidx;
  
  for(int j=0; j < causalcnt; j++) {
    config_new[j+1] = config[j];
  }
  
  config_new.sort();
  
  return config_new;
}


// swap a causal index snpidx with causal index j in config
IntegerVector swap(IntegerVector config, int snpidx, int j, int causalcnt) {
  
  IntegerVector config_new(causalcnt);
  
  for(int j=0; j < causalcnt; j++) {
    config_new[j] = config[j];
  }
  
  config_new[j] = snpidx;
  
  config_new.sort();
  
  return config_new;
}

// remove a causal index snpidx from config
IntegerVector remove(IntegerVector config, int snpidx, int causalcnt) {
  
  IntegerVector config_new(causalcnt-1);
  
  bool removed = false;
  
  for(int j=0; j < causalcnt; j++) {
    
    if(snpidx != config[j]) {
      
      if(!removed) {
        config_new[j] = config[j]; 
      } else {
        config_new[j-1] = config[j];
      }
      
    } else {
      removed = true;
    }
    
  }
  
  return config_new;
}


// generate neighborhood of the current config
// [[Rcpp::export]]
List nbh(IntegerVector config, int m) 
{
  List perm;
  
  int causalcnt = config.size();
  
  if(causalcnt==0 || m==1) {
    
    for(int i=0; i < m; i++) {
      
      perm.push_back(IntegerVector::create(i));
    }
    
  } else if(causalcnt==1) {
    
    int causalIdx = config[0];
    
    for(int i=0; i < m; i++) {
      
      if(i != causalIdx) {
        
        // swap
        perm.push_back(IntegerVector::create(i));
        
        // add
        IntegerVector config_new(2);
        
        if(i > causalIdx) {
          config_new(0) = causalIdx;
          config_new(1) = i;
        } else {
          config_new(0) = i;
          config_new(1) = causalIdx;
        }
        
        perm.push_back(config_new);
      }
    }
    
  } else if(causalcnt < m) {
    
    for(int i=0; i < m; i++) {
      
      bool hasIt = false;
      
      for(int j=0; j < causalcnt; j++) {
        
        if(i == config[j]) {
          hasIt = true;
        }
      }
      
      if(!hasIt) {
        
        // add
        perm.push_back(add(config, i, causalcnt));
        
        // swap
        for(int j=0; j < causalcnt; j++) {
          
          perm.push_back(swap(config, i, j, causalcnt));
        }
      } else {
        
        perm.push_back(remove(config, i, causalcnt));
      }
      
    }
    
  } else if (causalcnt == m && m > 1){ // causalcnt == m
    
    for(int i=0; i < m; i++) {
      perm.push_back(remove(config, i, causalcnt));
    }
  } else {
    throw std::runtime_error("nbh: unknown error");
  }
  
  return perm;
}

















