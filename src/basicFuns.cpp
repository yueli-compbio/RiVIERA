#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
vec dbetaPT(vec x, double alpha)
{
  return alpha * pow(x, alpha - 1);
}

vec lbetaPT(vec x, double alpha)
{
  return trunc_log(alpha) + (alpha - 1) * trunc_log(x);
}

double sigmoid(double x) {
  return 1/(1+trunc_exp(-x));
}

mat sigmoid(const mat& x) {
  return 1/(1+trunc_exp(-x));
}

// [[Rcpp::export]]
vec inferPri(const mat& ann, const vec& ann_w, double ann_w0)
{
  return sigmoid(ann * ann_w + ann_w0);
}

// [[Rcpp::export]]
vec inferPriByLoci(const mat& ann, 
                   const vec& ann_w, 
                   const vec& ann_w0, 
                   const umat& locusCursor)
{
  int global_up = 0;
  int global_dw = 0;
  
  vec pri(ann.n_rows);
  
  for(int l=0; l<(int)locusCursor.n_rows; l++) {
    
    global_up = locusCursor(l, 0);
    
    global_dw = locusCursor(l, 1);
    
    pri(span(global_up, global_dw)) = 
      
      inferPri(ann.rows(span(global_up, global_dw)), ann_w, ann_w0(l));
  }
  
  return pri;
}


double mod(int x, int m) {
  return x - m*floor(x/m);
}

// [[Rcpp::export]]
mat rwishart(double nu,  mat V) {
  // function to draw from Wishart (nu,V) and IW
  //
  // W ~ W(nu,V) E[W]=nuV
  //
  // WI=W^-1 E[WI]=V^-1/(nu-m-1)
  
  //    RNGScope rngscope;
  
  int m = V.n_rows;
  
  if (nu <= m - 1) {
    throw std::range_error("Inadmissible nu");
  }
  
  // Can't vectorise because Rcpp-sugar rchisq doesnt 
  // vectorise df argument
  // T has sqrt chisqs on diagonal and normals 
  // below diagonal and garbage above diagonal
  mat Z = randn<mat>(m, m);
  
  for(int i = 0; i < m; i++) {
    Z(i,i) = sqrt(R::rchisq(nu - i));
  }
  
  // explicitly declare T as triangular
  // top triangular is NaN ###
  mat C = trimatl(Z).t() * chol(V);
  
  // C is the upper triangular root of Wishart
  // therefore, W=C'C  this is the LU decomposition
  
  return C.t() * C;
}

// [[Rcpp::export]]
bool isPosDef(const mat& x) {
  
  vec eigval(x.n_rows);
  
  try{
    
    eigval = eig_sym(x);
    
  } catch(...) {
    
    // Rcout << "x: " << endl << x << endl;
    // throw std::runtime_error("check x above");
    
    return false;
  }
  
  return all(eigval > 1e-8);
}


// double dmvnrm_pinv(const vec& x, double mu, 
// mat sigma, bool logd = true)
// {
//   
//   int d = x.n_elem;
//   
//   return as_scalar(
//     -0.5 * d * log2pi - 0.5 * log(det(sigma)) -
//       0.5 * trans(x - mu) * pinv(sigma) * (x - mu));
// }


// [[Rcpp::export]]
double dmvnrm(const vec& x, double mu, 
              mat sigma, bool logd = true)
{
  int d = x.n_elem;
  
  mat rooti;
  
  // check positive-definiteness
  if(isPosDef(sigma)) {  
    
    try {
      
      rooti = trans(inv(trimatu(chol(sigma))));
      
    } catch(...) {
      
      return -INFINITY;
    }
    
  } else {
    
    return -INFINITY;
  }
    
  double rootisum = sum(log(rooti.diag()));
  
  double constants = -(static_cast<double>(d)/2.0) * log2pi;
  
  vec z = rooti * ( x - mu );
  
  double out = constants - 0.5 * sum(z % z) + rootisum;
  
  if (logd == false) {
    out = trunc_exp(out);
  }
  
  if(!is_finite(out)) {
    return -INFINITY;
  }
  
  return out;
}

// [[Rcpp::export]]
double dmvnrm_lambda(const vec& x, double mu, 
                     mat lambda, bool logd = true)
{
  return dmvnrm(x, mu, inv_sympd(lambda), logd);
}



// [[Rcpp::export]]
mat fastInv(mat x) {
  
  mat rooti = inv(trimatu(chol(x)));
  
  return rooti * rooti.t();
}


// [[Rcpp::export]]
vec duvnrm(const vec& x, double mu=0, double s2=1)
{
  return (1/sqrt((2.0*M_PI*s2))) * exp(-(1/(2*s2))*square(x-mu));
}

// [[Rcpp::export]]
umat createLocusCursor(const List& ldmat_list, int totalSNPs)
{
  int locNum = ldmat_list.size();
  
  umat locusCursor(locNum, 2);
  
  int cursor = 0;
  int m = 0;
  
  for(int l=0; l<locNum; l++) {
    
    m = as<mat>(ldmat_list[l]).n_rows;
    
    locusCursor(l,0) = cursor;
    locusCursor(l,1) = cursor + m - 1;
    
    cursor += m;
  }
  
  if(totalSNPs - 1 != (int)locusCursor(locusCursor.n_rows - 1, 1)) {
    
    Rcout << "totalSNPs: " << totalSNPs << endl;
    
    Rcout << "locusCursor(locusCursor.n_rows - 1, 1): ";
    Rcout << locusCursor(locusCursor.n_rows - 1, 1) << endl;
    
    throw std::invalid_argument(
        "Inconsistent number of SNPs in LD and other inputs");
  }
  
  return locusCursor;
}


// [[Rcpp::export]]
vec normalizeZscoreByLoci(const vec& zscore, const umat& locusCursor)
{
  int global_up = 0;
  
  int global_dw = 0;
  
  vec znorm(zscore.n_elem);
  
  for(int l=0; l<(int)locusCursor.n_rows; l++) {
    
    global_up = locusCursor(l, 0);
    
    global_dw = locusCursor(l, 1);
    
    vec zscore2 = square(zscore(span(global_up, global_dw)));
    
    znorm(span(global_up, global_dw)) = zscore2/accu(zscore2);
  }
  
  return znorm;
}






















