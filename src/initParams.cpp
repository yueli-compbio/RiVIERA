#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List initParams(const List& zscoreList, 
                const List& annList, 
                const List& ldmatList, 
                const List& locusCursorList,
                double causalCntPerLocusPrior,
                double s2_init,
                double pri0=1e-2) 
{
  int D = zscoreList.size();
  
  // initialize pip for each locus
  List pipList(D);
  
  // initialize prior for each locus as 1/locusSize
  List priList(D);
  
  List ann_w0(D);
  
  List ann_w0_tau(D);
  
  List ann_w0_mu(D);
  
  List s2(D);
  
  int K = 0;
  
  for(int d=0; d < D; d++) {
    
    vec zscore_d = zscoreList[d];
    
    vec pip_d = zeros<vec>(zscore_d.n_elem);
    
    mat ann_d = annList[d];
    
    if(K == 0) {
      
      K = ann_d.n_cols;
      
    } else if(K != (int)ann_d.n_cols) {
      throw std::invalid_argument("Inconsistent number of annotations");
    }
    
    mat locusCursor_d = locusCursorList[d];
    
    List ldmat_d = ldmatList[d];
    
    int L_d = ldmat_d.size();
    
    // cursor location of start and end snps per locus
    int global_up = 0;
    int global_dw = 0;
    
    // total SNPs for trait d
    int M_d = 0;
    
    vec pri_d(L_d);
    
    vec s2_d(L_d);
    
    if(ann_d.n_rows - 1 != locusCursor_d(locusCursor_d.n_rows - 1, 1)) {
      
      Rcout << "ann_d.n_rows: " << ann_d.n_rows << endl;
      
      Rcout << "locusCursor_d(locusCursor_d.n_rows - 1, 1): ";
      Rcout << locusCursor_d(locusCursor_d.n_rows - 1, 1) << endl;
      
      throw std::invalid_argument(
          "Last locusCursor must be last row of ann");
    }
    
    
    for(int l=0; l < L_d; l++) {
      
      global_up = locusCursor_d(l, 0);
      global_dw = locusCursor_d(l, 1);
      
      double locusSize_ld = global_dw - global_up + 1;
      
      mat ldmat_ld = ldmat_d[l];
      
      if(ldmat_ld.n_rows != locusSize_ld) {
        
        Rcout << "ldmat_ld.n_rows = " << ldmat_ld.n_rows << endl;
        
        Rcout << "locusSize_ld = " << locusSize_ld << endl;
        
        throw std::invalid_argument("Incorrect locusCursor");
      }
      
      M_d += locusSize_ld;
      
      // prevent log(0) below when there is only one snp in the locus
      if(locusSize_ld>1) {
        pri_d(l) = causalCntPerLocusPrior/locusSize_ld;  
      } else {
        pri_d(l) = pri0;
      }
      
      pip_d(span(global_up,global_dw)).fill(pri_d(l));
    }
    
    pri_d.elem(find(pri_d > pri0)).fill(pri0);
    
    priList[d] = pri_d;
    
    pipList[d] = pip_d;
    
    ann_w0_mu[d] = log(pri_d) - log(1-pri_d);
    
    ann_w0_tau[d] = zeros<vec>(L_d);
    
    ann_w0[d] = ann_w0_mu[d];
    
    if((int)zscore_d.n_elem != M_d) {
      throw std::invalid_argument("Incorrect input dimensions");
    }
    
    
    // initialize s2 the same for all loci all traits
    s2_d.fill(s2_init);
    
    s2[d] = s2_d;
  }
  
  // initialize model parameters
  mat ann_w = 1e-8*randn<mat>(K, D);
  
  mat ann_w_lambda = eye<mat>(D,D); // precision matrix for w_ik across D
  
  mat ann_w_sigma0 = eye(D, D);
  
  return List::create(Named("pipList")=pipList, 
                      Named("priList")=priList,
                      
                      Named("ann_w")=ann_w,
                      Named("ann_w_lambda")=ann_w_lambda,
                      Named("ann_w_sigma0")=ann_w_sigma0,
                      
                      Named("ann_w0")=ann_w0,
                      Named("ann_w0_mu")=ann_w0_mu,
                      Named("ann_w0_tau")=ann_w0_tau,
                      
                      Named("s2")=s2);
}






























