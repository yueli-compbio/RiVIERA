#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "calGrad.h"
#include "lpr.h"
#include "basicFuns.h"

double kineticsW(vec p_ann_w, vec p_ann_w0)
{
  return (accu(square(p_ann_w)) + accu(square(p_ann_w0)))/2;
}

// [[Rcpp::export]]
vec hmcW(const mat& ann,
          const List& configs_list, 
          const List& configs_post_list,
          const umat& locusCursor,
          vec ann_w_cur,
          vec ann_w_mu,
          mat ann_w_lambda,
          vec ann_w0_cur,
          vec ann_w0_mu,
          double ann_w0_tau,
          double step, int nsteps, 
          bool verbose)
{
  // save the current state in case rejection of the prposed state
  vec ann_w = ann_w_cur;
  
  vec ann_w0 = ann_w0_cur;
  
  int locNum = locusCursor.n_rows;
  
  int K = ann.n_cols;
  
  // initialize momentum variables
  vec p_ann_w = randn<vec>(K);
  vec p_ann_w0 = randn<vec>(locNum);
  
  // vec p_ann_w = zeros<vec>(K);
  // vec p_ann_w0 = zeros<vec>(locNum);
  
  double kinetics_cur = kineticsW(p_ann_w, p_ann_w0);
  
  vec grad = calGradW(
    
    ann, 
    
    configs_list,
    configs_post_list,
    locusCursor,
    
    ann_w, 
    ann_w_mu,
    ann_w_lambda, 
    
    ann_w0,
    ann_w0_mu,
    ann_w0_tau);
  
  // set momentum variables
  p_ann_w0 = p_ann_w0 - step * grad(span(0,locNum-1))/2;
  
  p_ann_w = p_ann_w - step * grad(span(locNum, locNum+K-1))/2;
  
  
  // Alternate full steps for position and momentum
  for (int iter=0; iter < nsteps; iter++) {
    
    // update position variables
    ann_w = ann_w + step * p_ann_w;
    
    ann_w0 = ann_w0 + step * p_ann_w0;
    
    // update gradients
    grad = calGradW(
      
      ann, 
      
      configs_list,
      configs_post_list,
      locusCursor,
      
      ann_w, 
      ann_w_mu,
      ann_w_lambda, 
      
      ann_w0,
      ann_w0_mu,
      ann_w0_tau);
    
    // update momentum variables with one full step
    p_ann_w0 = p_ann_w0 - step * grad(span(0,locNum-1));
    
    p_ann_w = p_ann_w - step * grad(span(locNum, locNum+K-1));
  }
  
  // approximate the final momentum
  grad = calGradW(
    
    ann, 
    
    configs_list,
    configs_post_list,
    locusCursor,
    
    ann_w, 
    ann_w_mu,
    ann_w_lambda, 
    
    ann_w0,
    ann_w0_mu,
    ann_w0_tau);
  
  // complete momentum variables with another half step
  p_ann_w0 = p_ann_w0 - step * grad(span(0,locNum-1))/2;
  
  p_ann_w = p_ann_w - step * grad(span(locNum, locNum+K-1))/2;
  
  
  // Evaluate potential and kinetic energies at the end of trajectory
  // hamiltonian dynamic of the proposed state
  double potential_new = -lprW(
    
    ann,
    
    configs_list,
    configs_post_list,
    locusCursor,
    
    ann_w, 
    ann_w_mu,
    ann_w_lambda,
    
    ann_w0,
    ann_w0_mu,
    ann_w0_tau);
  
  
  double kinetics_new = kineticsW(p_ann_w, p_ann_w0);
  
  double h_new = potential_new + kinetics_new;
  
  double potential_cur = -lprW(
    
    ann,
    
    configs_list,
    configs_post_list,
    locusCursor,
    
    ann_w_cur, 
    ann_w_mu,
    ann_w_lambda,
    
    ann_w0_cur,
    ann_w0_mu,
    ann_w0_tau);
  
  // hamiltonian dynamic of the current state
  double h_cur = potential_cur + kinetics_cur;
  
  // exp(-(H_new - H_cur))
  double dH = trunc_exp(h_cur - h_new);
  
  vec ret(locNum + K);
  
  if(randu() < dH) { // accept
    
    ret(span(0,locNum-1)) = ann_w0;
    
    ret(span(locNum, locNum+K-1)) = ann_w;
    
  } else {
    
    ret(span(0,locNum-1)) = ann_w0_cur;
    
    ret(span(locNum, locNum+K-1)) = ann_w_cur;
  }
  
  if(verbose) {
    Rcout << endl << "hmcW: " << endl;
    Rcout << "dH: " << dH << endl;
    Rcout << "potential_cur: " << potential_cur << endl;
    Rcout << "potential_new: " << potential_new << endl;
    Rcout << "kinetics_cur: " << kinetics_cur << endl;
    Rcout << "kinetics_new: " << kinetics_new << endl << endl;
    
    Rcout << "ann_w: " << ann_w << endl;
    Rcout << "ann_w0: " << ann_w0 << endl;
  }
  
  return ret;
}







