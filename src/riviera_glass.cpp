#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "basicFuns.h"
#include "riviera_ridge.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec calGradPT_glass(vec pri,
                    vec ppa,
                    mat ann,
                    vec ann_w,
                    double ann_w0,
                    double ann_w0_mu,
                    mat groupMat,
                    vec sdf)
{
  int K = ann.n_cols;
  
  vec grad(K+1);
  
  // intercept
  grad(0) = accu(ppa - pri) - (ann_w0 - ann_w0_mu);
  
  // weights
  grad(span(1,K)) = ann.t() * (ppa - pri) - 
    
    ((groupMat.t() * sdf) % ann_w) / (groupMat.t() * sqrt(groupMat * square(ann_w)));
  
  return grad;
}


// [[Rcpp::export]]
mat calHessPT_glass(vec pri, mat ann, vec ann_w, mat groupMat, vec sdf)
{
  int K = ann.n_cols;
  
  mat hess = zeros<mat>(K+1,K+1);
  
  hess(0,0) = as_scalar(pri.t() * (pri - 1)) - 1; // intercept only
  
  // intercept
  hess(span(1,K), span(0,0)) = ann.t() * (pri % (pri - 1));
  
  hess(span(0,0), span(1,K)) = hess(span(1,K), span(0,0)).t();
  
  // repeat K_g for ann in the same group g
  // repmat(KxG (GxK Kx1), 1, K) -> KxK
  mat ann_w_groupNorm2 = repmat(groupMat.t() * (groupMat * square(ann_w)), 1, K);
  
  mat groupenalty = repmat(groupMat.t() * sdf, 1, K) % 
  
  (diagmat(ones(K)) - (ann_w * ann_w.t())/ann_w_groupNorm2)/sqrt(ann_w_groupNorm2);
  
  // only annotations in the group count
  groupenalty = (groupMat.t() * groupMat) % groupenalty;
  
  hess(span(1,K),span(1,K)) = ann.t() * (ann % repmat(pri % (pri - 1), 1, K)) - groupenalty;
  
  return hess;
}


// [[Rcpp::export]]
double lprPT_glass(vec pval,
                   vec pri,
                   mat ann,
                   double alpha,
                   vec ann_w,
                   double ann_w0,
                   double ann_w0_mu,
                   mat groupMat,
                   vec sdf)
{
  // vec pri = inferPri(ann, ann_w, ann_w0);
  
  return accu((trunc_log(pri) + lbetaPT(pval, alpha)) + trunc_log(1-pri)) -
    accu(sdf % sqrt(groupMat * square(ann_w))) - 
    pow(ann_w0 - ann_w0_mu, 2)/2;
}


// [[Rcpp::export]]
double lprPT_glass_givenPPA(vec pval,
                            vec pri,
                            vec ppa,
                            double alpha,
                            vec ann_w,
                            double ann_w0,
                            double ann_w0_mu,
                            mat groupMat,
                            vec sdf)
{
  return accu( ppa % (trunc_log(pri) + lbetaPT(pval, alpha)) + 
               (1-ppa) % trunc_log(1-pri) ) - 
               accu(sdf % sqrt(groupMat * square(ann_w))) - 
               pow(ann_w0 - ann_w0_mu, 2)/2;
}


// [[Rcpp::export]]
void printAnn(std::vector< std::string > x, uvec activeGroupIdx, vec groupScores) {
  
  uvec q = sort_index(groupScores, "acscent");
  
  uvec isActive = zeros<uvec>(groupScores.n_elem);
  
  isActive(activeGroupIdx).ones();
  
  Rcout << endl;
  
  for( int i=0; i < (int)x.size(); i++ ) {
    
    int g = q(i);
    
    Rcout << x[g] << ": S(g) = " << groupScores(g) 
      << "; active = " << isActive(g) << endl;
  }
  
  Rcout << endl;
}


// [[Rcpp::export]]
List riviera_glass(vec pval,
                  const mat& ann,
                  mat groupMat,
                  uvec targetGroups,
                  std::vector<std::string> targetGroupNames,
                  double minPval=1e-30,
                  double alpha0=0.1,
                  double pri0 = 1e-3,
                  double epsilon=1e-2,
                  double thres=1e-3,
                  int max_iter=1e3,
                  double minAlpha=1e-30,
                  double maxAlpha=0.5,
                  double groupScoreThres=0.1,
                  int verbose=1)
{
  double alpha = alpha0;
  
  List emfit;
  
  vec lpr_init;
  
  uvec smallPval_idx = find(pval < minPval);
  
  int smallPval_cnt = smallPval_idx.n_elem;
  
  if(smallPval_cnt > 0) {
    
    Rcout << smallPval_cnt << 
      " SNPs have p-value < " << minPval << endl;
    
    Rcout << "They are set to " << minPval;
    Rcout << " for numerical stability" << endl;
    
    pval(smallPval_idx).fill(minPval);  
  }
  
  // Full EM update
  Rcout << "Estimate parameters by model by EM:" << endl;
  
  // intialize
  int N = pval.n_elem;
  int G = groupMat.n_rows;
  int K = ann.n_cols;
  
  double ann_w0_mu = trunc_log(pri0) - trunc_log(1-pri0);
  
  double ann_w0 = ann_w0_mu;
  
  vec ann_w = zeros<vec>(K);
  
  int iter = 0;
  
  double lpr_diff = 1;
  
  vec lpr = zeros<vec>(max_iter);
  
  // E-step
  vec pri = inferPri(ann, ann_w, ann_w0);
  
  vec causal = pri % dbetaPT(pval, alpha);
  
  vec ppa = causal/(causal + (1-pri));
  
  double N1 = accu(ppa); // causal SNPs count
  
  vec sdf = ones<vec>(G);
  
  // Evaluate
  lpr(iter) = lprPT_glass(pval,pri,ann,alpha,ann_w,ann_w0,ann_w0_mu,groupMat,sdf);
  
  Rcout << "iter " << iter+1 << ": lpr = " << lpr(iter);
  Rcout << "; K = " << K << "; G = " << G;
  Rcout << "; alpha = " << alpha << "; pri0 = " << 1/(1+exp(-ann_w0));
  Rcout << "; diff = " << lpr_diff << endl;
  
  iter++;
  
  vec lpval = trunc_log(pval);
  
  mat ann_grp = ann * groupMat.t();
  
  ann_grp(find(ann_grp>0)).ones();
  
  vec ann_w_grp = 1e-5 * randn<vec>(ann_grp.n_cols);
  
  mat groupMat_tar = groupMat.rows(targetGroups);
  
  vec targetGroupScores = ann_w_grp(targetGroups);
  
  int targetGroupNum = targetGroups.n_elem;
  
  vec ann_w_mu = zeros<vec>(G);
  
  while(iter < max_iter && (iter < 10 || lpr_diff > thres)) {
    
    // M-step
    alpha = -N1/accu(ppa % lpval);
    
    if(alpha > maxAlpha) {
      alpha = maxAlpha;
    } else if (alpha < minAlpha) {
      alpha = minAlpha;
    }
    
    mat hess_grp = calHessPT_ridge(pri, ann_grp);
    
    vec newtonUpdates = pinv(hess_grp) * 
      calGradPT_ridge(pri, ppa, ann_grp, ann_w_grp, ann_w_mu, ann_w0, ann_w0_mu);
    
    ann_w0 -= epsilon * newtonUpdates(0);
    
    ann_w_grp -= epsilon * newtonUpdates(span(1,G));

    // covariance is the inverse of negative hessian matrix
    // target group z-score
    mat ann_w_cov = pinv(-hess_grp);
    
    vec ann_w_se = sqrt(abs(ann_w_cov.diag()));
    
    targetGroupScores = ann_w_grp(targetGroups)/ann_w_se(targetGroups);
    
    // targetGroupScores = square(targetGroupScores);
    
    uvec activeGroupIdx = find(targetGroupScores/max(targetGroupScores) > groupScoreThres);
    
    uvec inActiveGroupIdx = find(targetGroupScores/max(targetGroupScores) <= groupScoreThres);
    
    if(verbose==2) {
      printAnn(targetGroupNames, activeGroupIdx, targetGroupScores);  
    }
    
    ann_w(span(0, G-targetGroupNum)) = ann_w_grp(span(0, G-targetGroupNum));
    
    pri = inferPri(ann, ann_w, ann_w0);
    
    for(int j=0; j<(int)activeGroupIdx.n_elem; j++) {
      
      int g = activeGroupIdx(j);
      
      uvec annIdx_g = find(groupMat_tar.row(g)==1);
      
      vec ann_w_g = ann_w(annIdx_g);
      
      int K_g = ann_w_g.n_elem;
      
      vec grad_w = ann.cols(annIdx_g).t() * (ppa - pri);
      
      // add a small step to ensure differentiable weights
      if(accu(ann_w_g)==0) {
        ann_w_g += 1e-8 * grad_w;
      }
      
      grad_w -= sdf(g) * ann_w_g / sqrt(accu(square(ann_w_g)));
      
      double ann_w_g_norm2 = accu(square(ann_w_g));
      
      mat hess_w = ann.cols(annIdx_g).t() * (ann.cols(annIdx_g) % repmat(pri % (pri - 1), 1, K_g)) -
        sdf(g) * (diagmat(ones(K_g)) - (ann_w_g * ann_w_g.t())/ann_w_g_norm2)/sqrt(ann_w_g_norm2);
      
      ann_w(annIdx_g) -= epsilon * (pinv(hess_w) * grad_w);
    }
    
    for(int g=0; g<(int)inActiveGroupIdx.n_elem; g++) {
      
      ann_w(find(groupMat_tar.row(inActiveGroupIdx(g))==1)).zeros();
    }
    
    // Evaluate
    lpr(iter) = lprPT_glass(pval, inferPri(ann, ann_w, ann_w0),
        ann, alpha, ann_w, ann_w0, ann_w0_mu, groupMat, sdf);
    
    // lpr(iter) = lprPT_ridge(pval, pri, ann_grp, alpha, 
    //     ann_w_grp, ann_w_mu, ann_w0, ann_w0_mu);
    
    lpr_diff = (lpr(iter) - lpr(iter-1))/lpr(iter);
    
    if(lpr_diff<0) {
      lpr_diff = -lpr_diff;
    }
    
    Rcout << "iter " << iter+1 << ": lpr = " << lpr(iter);
    Rcout << "; alpha = " << alpha << "; pri0 = " << 1/(1+exp(-ann_w0));
    Rcout << "; diff = " << lpr_diff << endl;
    
    // E-step
    pri = inferPri(ann_grp, ann_w_grp, ann_w0);
    
    causal = pri % dbetaPT(pval, alpha);
    
    ppa = causal/(causal + (1-pri));
    
    N1 = accu(ppa); // causal SNPs count
    
    iter++;
  }
  
  lpr = lpr(span(0,iter-1));
  
  
  // group-level z-score
  // covariance is the inverse of negative hessian matrix
  mat ann_w_cov_grp = pinv(-calHessPT_ridge(inferPri(ann_grp, ann_w_grp, ann_w0), ann_grp));
  vec ann_w_se_grp = sqrt(abs(ann_w_cov_grp.diag()));
  vec ann_w_zscore_grp = ann_w_grp/ann_w_se_grp(span(1,G));
  
  
  // target annotations z-score
  // covariance is the inverse of negative hessian matrix
  mat ann_w_cov = diagmat(ones<vec>(K));
  
  pri = inferPri(ann, ann_w, ann_w0);
  
  uvec isActiveAnnot = zeros<uvec>(K);
    
  isActiveAnnot(find(ann_w!=0)).ones();
  
  uvec activeAnnotIdx = find(isActiveAnnot==1);
  
  int K_active = activeAnnotIdx.n_elem;
  
  uvec activeGroupIdx = find(sum(groupMat.cols(activeAnnotIdx),1) > 0);
  
  mat inv_nhess = 
    pinv(-calHessPT_glass(pri, ann.cols(activeAnnotIdx), ann_w(activeAnnotIdx), 
                          groupMat(activeGroupIdx,activeAnnotIdx), sdf(activeGroupIdx)));
  
  ann_w_cov(activeAnnotIdx, activeAnnotIdx) = inv_nhess(span(1,K_active), span(1,K_active));
  
  vec ann_w_se = sqrt(abs(ann_w_cov.diag()));
  
  vec ann_w_zscore = ann_w/ann_w_se;
  
  
  return List::create(Named("alpha")=alpha,
                      Named("pri0")=pri0,
                      Named("ann_w_grp")=ann_w_grp,
                      Named("ann_w_cov_grp")=ann_w_cov_grp,
                      Named("ann_w_se_grp")=ann_w_se_grp,
                      Named("ann_w_zscore_grp")=ann_w_zscore_grp,
                      Named("ann_w")=ann_w,
                      Named("ann_w0")=ann_w0,
                      Named("ann_w_se")=ann_w_se,
                      Named("ann_w_cov")=ann_w_cov,
                      Named("ann_w_zscore")=ann_w_zscore,
                      Named("lpr")=lpr);
}
























