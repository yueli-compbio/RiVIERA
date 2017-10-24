#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double lprW_l(const vec& pri,
              const umat& configs,
              const vec& configs_post,
              
              const vec& ann_w,
              const vec& ann_w_mu,
              const mat& ann_w_lambda,
              
              double ann_w0,
              double ann_w0_mu,
              double ann_w0_tau);

double lprW(const mat& ann,
            const List& configs_list,
            const List& configs_post_list,
            
            const umat& locusCursor,
            
            const vec& ann_w,
            const vec& ann_w_mu,
            const mat& ann_w_lambda,
            
            vec ann_w0,
            vec ann_w0_mu,
            double ann_w0_tau);

double lprS_l(const vec& zscore,
              const mat& ldmat,
              const umat& configs,
              const vec configs_post,
              double s2);

double lprS(const vec& zscore,
            const List& ldmat_list,
            const List& configs_list,
            const List& configs_post_list,
            const umat& locusCursor,
            const vec& s2);
