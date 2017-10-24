#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "calGrad.h"
#include "lpr.h"


// check and correct for both variable q and momentum p
// if q is outside of defined bounds
void check_bound(vec& q, vec& p, 
                 double lb, double ub, const char *qname) 
{
  
  bool outbound = true;
  
  // constrain lbBetaMeanVal < mu < ubBetaMeanVal
  while (outbound) {
    
    for (int i=0; i<(int)q.n_elem; i++) {
      
      outbound = true;
      
      if(q(i) < lb) {
        
        Rf_warning("%s = %.3f < lower bound %.3f!", qname, q(i), lb);
        
        q(i) = lb + (lb - q(i));
        
      } else if (q(i) > ub) {
        
        Rf_warning("%s = %.3f > upper bound %.3f!", qname, q(i), ub);
        
        q(i) = ub - (q(i) - ub);
        
      } else {
        outbound = false;
      }
      
      if (outbound) {
        p(i) = -p(i);
      }
    }
  }
}

double kineticsS(const vec& p_s2) {
  
  return accu(square(p_s2))/2;
}

// [[Rcpp::export]]
vec hmcS(const vec& zscore, 
         const List& ldmat_list,
         const List& configs_list, 
         const List& configs_post_list,
         const umat& locusCursor,
         vec s2_cur,
         double step, int nsteps,
         bool verbose, 
         double lbs2 = 1e-8, 
         double ubs2 = 1e3)
{
  vec s2 = s2_cur;
  
  int locNum = s2.n_elem;
  
  vec p_s2 = randn<vec>(locNum);
  // vec p_s2 = zeros<vec>(locNum);
  
  double kinetics_cur = kineticsS(p_s2);
  
  // list to store gradients
  vec grad = calGradS(zscore,
                      ldmat_list,
                      configs_list,
                      configs_post_list,
                      locusCursor,
                      s2);
  
  p_s2 = p_s2 - step * grad/2;
  
  // Alternate full steps for position and momentum
  for (int iter=0; iter < nsteps; iter++) {
    
    s2 = s2 + step * p_s2;
    
    // update gradients
    grad = calGradS(zscore,
                    ldmat_list,
                    configs_list,
                    configs_post_list,
                    locusCursor,
                    s2);
    
    // update momentum variables with one full step
    p_s2 = p_s2 - step * grad;
    
    check_bound(s2, p_s2, lbs2, ubs2, "s2");
    
    // if(!is_finite(s2)) {
    //   
    //   Rcout << "configs_post_list: " << as<vec>(configs_post_list[0]) << endl;
    //   Rcout << "s2: " << s2 << endl;
    //   Rcout << "p_s2: " << p_s2 << endl;
    //   Rcout << "grad: " << grad << endl;
    //   
    //   throw std::runtime_error("s2 is not finite!");
    // }
  }
  
  // approximate the final momentum
  grad = calGradS(zscore,
                  ldmat_list,
                  configs_list,
                  configs_post_list,
                  locusCursor,
                  s2);
  
  // complete momentum variables with another half step
  p_s2 = p_s2 - step * grad/2;
  
  
  // Evaluate potential and kinetic energies at the end of trajectory
  
  // hamiltonian dynamic of the proposed state
  double potential_new = -lprS(zscore, 
                               ldmat_list, 
                               configs_list, 
                               configs_post_list, 
                               locusCursor, s2);
  
  double kinetics_new = kineticsS(p_s2);
  
  double h_new = potential_new + kinetics_new;
  
  double potential_cur = -lprS(zscore, 
                               ldmat_list, 
                               configs_list, 
                               configs_post_list, 
                               locusCursor, s2_cur);
  
  // hamiltonian dynamic of the current state
  double h_cur = potential_cur + kinetics_cur;
  
  // exp(-(H_new - H_cur))
  double dH = trunc_exp(h_cur - h_new);
  
  vec ret;
  
  if(randu() < dH) { // accept
    
    ret = s2;
    
  } else {
    
    ret = s2_cur;
  }
  

  if(verbose) {
    Rcout << endl << "hmcS: " << endl;
    Rcout << "dH: " << dH << endl;
    Rcout << "potential_cur: " << potential_cur << endl;
    Rcout << "potential_new: " << potential_new << endl;
    Rcout << "kinetics_cur: " << kinetics_cur << endl;
    Rcout << "kinetics_new: " << kinetics_new << endl;
    
    Rcout << "s2: " << s2 << endl;
  }
  
  if(!is_finite(dH)) {
    
    Rcout << "s2: " << s2 << endl;
    Rcout << "s2_cur: " << s2_cur << endl;
    Rcout << "configs_post_list: " << 
      as<vec>(configs_post_list[0]) << endl;
    Rcout << "configs_list: " << as<umat>(configs_list[0]) << endl;
    
    throw std::runtime_error("dH is not finite!");
  }
  
  
  
  return ret;
}











