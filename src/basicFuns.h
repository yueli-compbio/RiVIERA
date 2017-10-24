#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec dbetaPT(vec x, double alpha);

vec lbetaPT(vec x, double alpha);

double sigmoid(double x);

mat sigmoid(const mat& x);

vec inferPri(const mat& ann, const vec& ann_w, double ann_w0);

vec inferPriByLoci(const mat& ann, 
                   const vec& ann_w, 
                   const vec& ann_w0, 
                   const umat& locusCursor);

double mod(int x, int m);

mat rwishart(double nu,  mat V);

double dmvnrm(const vec& x, double mu, mat sigma, bool logd = true);

double dmvnrm_lambda(const vec& x, double mu, mat lambda, bool logd = true);

mat fastInv(mat x);

vec duvnrm(const vec& x, double mu=0, double s2=1);

umat createLocusCursor(const List& ldmat_list, int totalSNPs);

vec normalizeZscoreByLoci(const vec& zscore, const umat& locusCursor);



