#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// dbetaPT
vec dbetaPT(vec x, double alpha);
RcppExport SEXP _RiVIERA_dbetaPT(SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetaPT(x, alpha));
    return rcpp_result_gen;
END_RCPP
}
// inferPri
vec inferPri(const mat& ann, const vec& ann_w, double ann_w0);
RcppExport SEXP _RiVIERA_inferPri(SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    rcpp_result_gen = Rcpp::wrap(inferPri(ann, ann_w, ann_w0));
    return rcpp_result_gen;
END_RCPP
}
// inferPriByLoci
vec inferPriByLoci(const mat& ann, const vec& ann_w, const vec& ann_w0, const umat& locusCursor);
RcppExport SEXP _RiVIERA_inferPriByLoci(SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP locusCursorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    rcpp_result_gen = Rcpp::wrap(inferPriByLoci(ann, ann_w, ann_w0, locusCursor));
    return rcpp_result_gen;
END_RCPP
}
// rwishart
mat rwishart(double nu, mat V);
RcppExport SEXP _RiVIERA_rwishart(SEXP nuSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(rwishart(nu, V));
    return rcpp_result_gen;
END_RCPP
}
// isPosDef
bool isPosDef(const mat& x);
RcppExport SEXP _RiVIERA_isPosDef(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(isPosDef(x));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm
double dmvnrm(const vec& x, double mu, mat sigma, bool logd);
RcppExport SEXP _RiVIERA_dmvnrm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm(x, mu, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_lambda
double dmvnrm_lambda(const vec& x, double mu, mat lambda, bool logd);
RcppExport SEXP _RiVIERA_dmvnrm_lambda(SEXP xSEXP, SEXP muSEXP, SEXP lambdaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_lambda(x, mu, lambda, logd));
    return rcpp_result_gen;
END_RCPP
}
// fastInv
mat fastInv(mat x);
RcppExport SEXP _RiVIERA_fastInv(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fastInv(x));
    return rcpp_result_gen;
END_RCPP
}
// duvnrm
vec duvnrm(const vec& x, double mu, double s2);
RcppExport SEXP _RiVIERA_duvnrm(SEXP xSEXP, SEXP muSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(duvnrm(x, mu, s2));
    return rcpp_result_gen;
END_RCPP
}
// createLocusCursor
umat createLocusCursor(const List& ldmat_list, int totalSNPs);
RcppExport SEXP _RiVIERA_createLocusCursor(SEXP ldmat_listSEXP, SEXP totalSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type ldmat_list(ldmat_listSEXP);
    Rcpp::traits::input_parameter< int >::type totalSNPs(totalSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(createLocusCursor(ldmat_list, totalSNPs));
    return rcpp_result_gen;
END_RCPP
}
// normalizeZscoreByLoci
vec normalizeZscoreByLoci(const vec& zscore, const umat& locusCursor);
RcppExport SEXP _RiVIERA_normalizeZscoreByLoci(SEXP zscoreSEXP, SEXP locusCursorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    rcpp_result_gen = Rcpp::wrap(normalizeZscoreByLoci(zscore, locusCursor));
    return rcpp_result_gen;
END_RCPP
}
// estimateParByEM
List estimateParByEM(vec pval, double alpha_init, double pri0_init, double minPval, double thres, int max_iter);
RcppExport SEXP _RiVIERA_estimateParByEM(SEXP pvalSEXP, SEXP alpha_initSEXP, SEXP pri0_initSEXP, SEXP minPvalSEXP, SEXP thresSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_init(alpha_initSEXP);
    Rcpp::traits::input_parameter< double >::type pri0_init(pri0_initSEXP);
    Rcpp::traits::input_parameter< double >::type minPval(minPvalSEXP);
    Rcpp::traits::input_parameter< double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateParByEM(pval, alpha_init, pri0_init, minPval, thres, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// calGradW_l
vec calGradW_l(const mat& ann, const umat& configs, const vec& configs_post, const vec& ann_w, const vec& ann_w_mu, const mat& ann_w_lambda, double ann_w0, double ann_w0_mu, double ann_w0_tau);
RcppExport SEXP _RiVIERA_calGradW_l(SEXP annSEXP, SEXP configsSEXP, SEXP configs_postSEXP, SEXP ann_wSEXP, SEXP ann_w_muSEXP, SEXP ann_w_lambdaSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP ann_w0_tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const umat& >::type configs(configsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type configs_post(configs_postSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann_w_lambda(ann_w_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_tau(ann_w0_tauSEXP);
    rcpp_result_gen = Rcpp::wrap(calGradW_l(ann, configs, configs_post, ann_w, ann_w_mu, ann_w_lambda, ann_w0, ann_w0_mu, ann_w0_tau));
    return rcpp_result_gen;
END_RCPP
}
// calGradW
vec calGradW(const mat& ann, const List& configs_list, const List& configs_post_list, const umat& locusCursor, const vec& ann_w, const vec& ann_w_mu, const mat& ann_w_lambda, vec ann_w0, vec ann_w0_mu, double ann_w0_tau);
RcppExport SEXP _RiVIERA_calGradW(SEXP annSEXP, SEXP configs_listSEXP, SEXP configs_post_listSEXP, SEXP locusCursorSEXP, SEXP ann_wSEXP, SEXP ann_w_muSEXP, SEXP ann_w_lambdaSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP ann_w0_tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_list(configs_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_post_list(configs_post_listSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann_w_lambda(ann_w_lambdaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_tau(ann_w0_tauSEXP);
    rcpp_result_gen = Rcpp::wrap(calGradW(ann, configs_list, configs_post_list, locusCursor, ann_w, ann_w_mu, ann_w_lambda, ann_w0, ann_w0_mu, ann_w0_tau));
    return rcpp_result_gen;
END_RCPP
}
// calGrad_dEdV
mat calGrad_dEdV(vec x, mat sigma);
RcppExport SEXP _RiVIERA_calGrad_dEdV(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(calGrad_dEdV(x, sigma));
    return rcpp_result_gen;
END_RCPP
}
// calGrad_dEdS
double calGrad_dEdS(vec x, mat sigma, double s);
RcppExport SEXP _RiVIERA_calGrad_dEdS(SEXP xSEXP, SEXP sigmaSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(calGrad_dEdS(x, sigma, s));
    return rcpp_result_gen;
END_RCPP
}
// calGradS_l
double calGradS_l(const vec& zscore, const mat& ldmat, const umat& configs, const vec& configs_post, double s2);
RcppExport SEXP _RiVIERA_calGradS_l(SEXP zscoreSEXP, SEXP ldmatSEXP, SEXP configsSEXP, SEXP configs_postSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ldmat(ldmatSEXP);
    Rcpp::traits::input_parameter< const umat& >::type configs(configsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type configs_post(configs_postSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(calGradS_l(zscore, ldmat, configs, configs_post, s2));
    return rcpp_result_gen;
END_RCPP
}
// calGradS
vec calGradS(const vec& zscore, const List& ldmat_list, const List& configs_list, const List& configs_post_list, const umat& locusCursor, const vec& s2);
RcppExport SEXP _RiVIERA_calGradS(SEXP zscoreSEXP, SEXP ldmat_listSEXP, SEXP configs_listSEXP, SEXP configs_post_listSEXP, SEXP locusCursorSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const List& >::type ldmat_list(ldmat_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_list(configs_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_post_list(configs_post_listSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    Rcpp::traits::input_parameter< const vec& >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(calGradS(zscore, ldmat_list, configs_list, configs_post_list, locusCursor, s2));
    return rcpp_result_gen;
END_RCPP
}
// hmcS
vec hmcS(const vec& zscore, const List& ldmat_list, const List& configs_list, const List& configs_post_list, const umat& locusCursor, vec s2_cur, double step, int nsteps, bool verbose, double lbs2, double ubs2);
RcppExport SEXP _RiVIERA_hmcS(SEXP zscoreSEXP, SEXP ldmat_listSEXP, SEXP configs_listSEXP, SEXP configs_post_listSEXP, SEXP locusCursorSEXP, SEXP s2_curSEXP, SEXP stepSEXP, SEXP nstepsSEXP, SEXP verboseSEXP, SEXP lbs2SEXP, SEXP ubs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const List& >::type ldmat_list(ldmat_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_list(configs_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_post_list(configs_post_listSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    Rcpp::traits::input_parameter< vec >::type s2_cur(s2_curSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type lbs2(lbs2SEXP);
    Rcpp::traits::input_parameter< double >::type ubs2(ubs2SEXP);
    rcpp_result_gen = Rcpp::wrap(hmcS(zscore, ldmat_list, configs_list, configs_post_list, locusCursor, s2_cur, step, nsteps, verbose, lbs2, ubs2));
    return rcpp_result_gen;
END_RCPP
}
// hmcW
vec hmcW(const mat& ann, const List& configs_list, const List& configs_post_list, const umat& locusCursor, vec ann_w_cur, vec ann_w_mu, mat ann_w_lambda, vec ann_w0_cur, vec ann_w0_mu, double ann_w0_tau, double step, int nsteps, bool verbose);
RcppExport SEXP _RiVIERA_hmcW(SEXP annSEXP, SEXP configs_listSEXP, SEXP configs_post_listSEXP, SEXP locusCursorSEXP, SEXP ann_w_curSEXP, SEXP ann_w_muSEXP, SEXP ann_w_lambdaSEXP, SEXP ann_w0_curSEXP, SEXP ann_w0_muSEXP, SEXP ann_w0_tauSEXP, SEXP stepSEXP, SEXP nstepsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_list(configs_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_post_list(configs_post_listSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w_cur(ann_w_curSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< mat >::type ann_w_lambda(ann_w_lambdaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w0_cur(ann_w0_curSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_tau(ann_w0_tauSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(hmcW(ann, configs_list, configs_post_list, locusCursor, ann_w_cur, ann_w_mu, ann_w_lambda, ann_w0_cur, ann_w0_mu, ann_w0_tau, step, nsteps, verbose));
    return rcpp_result_gen;
END_RCPP
}
// inferConfigs_imp
List inferConfigs_imp(const vec& proposal, const vec& zscore, const mat& ldmat, const mat& ann, const vec& ann_w, double ann_w0, double s2, int sampleSize, bool verbose);
RcppExport SEXP _RiVIERA_inferConfigs_imp(SEXP proposalSEXP, SEXP zscoreSEXP, SEXP ldmatSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP s2SEXP, SEXP sampleSizeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ldmat(ldmatSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< int >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(inferConfigs_imp(proposal, zscore, ldmat, ann, ann_w, ann_w0, s2, sampleSize, verbose));
    return rcpp_result_gen;
END_RCPP
}
// inferConfigs_nbh
List inferConfigs_nbh(const vec& zscore, const mat& ldmat, const mat& ann, const vec& ann_w, double ann_w0, double s2, int sampleSize, bool verbose);
RcppExport SEXP _RiVIERA_inferConfigs_nbh(SEXP zscoreSEXP, SEXP ldmatSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP s2SEXP, SEXP sampleSizeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ldmat(ldmatSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< int >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(inferConfigs_nbh(zscore, ldmat, ann, ann_w, ann_w0, s2, sampleSize, verbose));
    return rcpp_result_gen;
END_RCPP
}
// inferConfigsAllLoci
List inferConfigsAllLoci(const vec& proposal, const vec& zscore, const List& ldmat_list, const umat& locusCursor, const mat& ann, const vec& ann_w, const vec& ann_w0, const vec& s2, int sampleSize, int sampleScheme, bool verbose);
RcppExport SEXP _RiVIERA_inferConfigsAllLoci(SEXP proposalSEXP, SEXP zscoreSEXP, SEXP ldmat_listSEXP, SEXP locusCursorSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP s2SEXP, SEXP sampleSizeSEXP, SEXP sampleSchemeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const List& >::type ldmat_list(ldmat_listSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< int >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< int >::type sampleScheme(sampleSchemeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(inferConfigsAllLoci(proposal, zscore, ldmat_list, locusCursor, ann, ann_w, ann_w0, s2, sampleSize, sampleScheme, verbose));
    return rcpp_result_gen;
END_RCPP
}
// inferPPA
vec inferPPA(vec pval, const mat& ann, const vec& ann_w, double ann_w0, double alpha, double minPval);
RcppExport SEXP _RiVIERA_inferPPA(SEXP pvalSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP alphaSEXP, SEXP minPvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type minPval(minPvalSEXP);
    rcpp_result_gen = Rcpp::wrap(inferPPA(pval, ann, ann_w, ann_w0, alpha, minPval));
    return rcpp_result_gen;
END_RCPP
}
// inferPPAByLoci
vec inferPPAByLoci(const vec& zscore, const mat& ann, const vec& ann_w, const vec& ann_w0, const vec& s2, const umat& locusCursor);
RcppExport SEXP _RiVIERA_inferPPAByLoci(SEXP zscoreSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP s2SEXP, SEXP locusCursorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< const vec& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    rcpp_result_gen = Rcpp::wrap(inferPPAByLoci(zscore, ann, ann_w, ann_w0, s2, locusCursor));
    return rcpp_result_gen;
END_RCPP
}
// initParams
List initParams(const List& zscoreList, const List& annList, const List& ldmatList, const List& locusCursorList, double causalCntPerLocusPrior, double s2_init, double pri0);
RcppExport SEXP _RiVIERA_initParams(SEXP zscoreListSEXP, SEXP annListSEXP, SEXP ldmatListSEXP, SEXP locusCursorListSEXP, SEXP causalCntPerLocusPriorSEXP, SEXP s2_initSEXP, SEXP pri0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type zscoreList(zscoreListSEXP);
    Rcpp::traits::input_parameter< const List& >::type annList(annListSEXP);
    Rcpp::traits::input_parameter< const List& >::type ldmatList(ldmatListSEXP);
    Rcpp::traits::input_parameter< const List& >::type locusCursorList(locusCursorListSEXP);
    Rcpp::traits::input_parameter< double >::type causalCntPerLocusPrior(causalCntPerLocusPriorSEXP);
    Rcpp::traits::input_parameter< double >::type s2_init(s2_initSEXP);
    Rcpp::traits::input_parameter< double >::type pri0(pri0SEXP);
    rcpp_result_gen = Rcpp::wrap(initParams(zscoreList, annList, ldmatList, locusCursorList, causalCntPerLocusPrior, s2_init, pri0));
    return rcpp_result_gen;
END_RCPP
}
// lprW
double lprW(const mat& ann, const List& configs_list, const List& configs_post_list, const umat& locusCursor, const vec& ann_w, const vec& ann_w_mu, const mat& ann_w_lambda, vec ann_w0, vec ann_w0_mu, double ann_w0_tau);
RcppExport SEXP _RiVIERA_lprW(SEXP annSEXP, SEXP configs_listSEXP, SEXP configs_post_listSEXP, SEXP locusCursorSEXP, SEXP ann_wSEXP, SEXP ann_w_muSEXP, SEXP ann_w_lambdaSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP ann_w0_tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_list(configs_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type configs_post_list(configs_post_listSEXP);
    Rcpp::traits::input_parameter< const umat& >::type locusCursor(locusCursorSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann_w_lambda(ann_w_lambdaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_tau(ann_w0_tauSEXP);
    rcpp_result_gen = Rcpp::wrap(lprW(ann, configs_list, configs_post_list, locusCursor, ann_w, ann_w_mu, ann_w_lambda, ann_w0, ann_w0_mu, ann_w0_tau));
    return rcpp_result_gen;
END_RCPP
}
// lprS_l
double lprS_l(const vec& zscore, const mat& ldmat, const umat& configs, const vec configs_post, double s2);
RcppExport SEXP _RiVIERA_lprS_l(SEXP zscoreSEXP, SEXP ldmatSEXP, SEXP configsSEXP, SEXP configs_postSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ldmat(ldmatSEXP);
    Rcpp::traits::input_parameter< const umat& >::type configs(configsSEXP);
    Rcpp::traits::input_parameter< const vec >::type configs_post(configs_postSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(lprS_l(zscore, ldmat, configs, configs_post, s2));
    return rcpp_result_gen;
END_RCPP
}
// nbh
List nbh(IntegerVector config, int m);
RcppExport SEXP _RiVIERA_nbh(SEXP configSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type config(configSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(nbh(config, m));
    return rcpp_result_gen;
END_RCPP
}
// riviera_fmap
List riviera_fmap(const vec& zscore, const List& ldmat_list, const mat& ann, const vec& ann_w_mu, const mat& ann_w_sigma0, int configSampleScheme, double c0, double pri0_init, double s2_init, double burninFrac, int sampleSize, int maxConfig, int max_iter, double step_W, int nsteps_W, double step_S, int nsteps_S, double thres, bool verbose, bool verboseC, bool verboseS, bool verboseW);
RcppExport SEXP _RiVIERA_riviera_fmap(SEXP zscoreSEXP, SEXP ldmat_listSEXP, SEXP annSEXP, SEXP ann_w_muSEXP, SEXP ann_w_sigma0SEXP, SEXP configSampleSchemeSEXP, SEXP c0SEXP, SEXP pri0_initSEXP, SEXP s2_initSEXP, SEXP burninFracSEXP, SEXP sampleSizeSEXP, SEXP maxConfigSEXP, SEXP max_iterSEXP, SEXP step_WSEXP, SEXP nsteps_WSEXP, SEXP step_SSEXP, SEXP nsteps_SSEXP, SEXP thresSEXP, SEXP verboseSEXP, SEXP verboseCSEXP, SEXP verboseSSEXP, SEXP verboseWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< const List& >::type ldmat_list(ldmat_listSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann_w_sigma0(ann_w_sigma0SEXP);
    Rcpp::traits::input_parameter< int >::type configSampleScheme(configSampleSchemeSEXP);
    Rcpp::traits::input_parameter< double >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< double >::type pri0_init(pri0_initSEXP);
    Rcpp::traits::input_parameter< double >::type s2_init(s2_initSEXP);
    Rcpp::traits::input_parameter< double >::type burninFrac(burninFracSEXP);
    Rcpp::traits::input_parameter< int >::type sampleSize(sampleSizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxConfig(maxConfigSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type step_W(step_WSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps_W(nsteps_WSEXP);
    Rcpp::traits::input_parameter< double >::type step_S(step_SSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps_S(nsteps_SSEXP);
    Rcpp::traits::input_parameter< double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseC(verboseCSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseS(verboseSSEXP);
    Rcpp::traits::input_parameter< bool >::type verboseW(verboseWSEXP);
    rcpp_result_gen = Rcpp::wrap(riviera_fmap(zscore, ldmat_list, ann, ann_w_mu, ann_w_sigma0, configSampleScheme, c0, pri0_init, s2_init, burninFrac, sampleSize, maxConfig, max_iter, step_W, nsteps_W, step_S, nsteps_S, thres, verbose, verboseC, verboseS, verboseW));
    return rcpp_result_gen;
END_RCPP
}
// calGradPT_glass
vec calGradPT_glass(vec pri, vec ppa, mat ann, vec ann_w, double ann_w0, double ann_w0_mu, mat groupMat, vec sdf);
RcppExport SEXP _RiVIERA_calGradPT_glass(SEXP priSEXP, SEXP ppaSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP groupMatSEXP, SEXP sdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pri(priSEXP);
    Rcpp::traits::input_parameter< vec >::type ppa(ppaSEXP);
    Rcpp::traits::input_parameter< mat >::type ann(annSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< mat >::type groupMat(groupMatSEXP);
    Rcpp::traits::input_parameter< vec >::type sdf(sdfSEXP);
    rcpp_result_gen = Rcpp::wrap(calGradPT_glass(pri, ppa, ann, ann_w, ann_w0, ann_w0_mu, groupMat, sdf));
    return rcpp_result_gen;
END_RCPP
}
// calHessPT_glass
mat calHessPT_glass(vec pri, mat ann, vec ann_w, mat groupMat, vec sdf);
RcppExport SEXP _RiVIERA_calHessPT_glass(SEXP priSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP groupMatSEXP, SEXP sdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pri(priSEXP);
    Rcpp::traits::input_parameter< mat >::type ann(annSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< mat >::type groupMat(groupMatSEXP);
    Rcpp::traits::input_parameter< vec >::type sdf(sdfSEXP);
    rcpp_result_gen = Rcpp::wrap(calHessPT_glass(pri, ann, ann_w, groupMat, sdf));
    return rcpp_result_gen;
END_RCPP
}
// lprPT_glass
double lprPT_glass(vec pval, vec pri, mat ann, double alpha, vec ann_w, double ann_w0, double ann_w0_mu, mat groupMat, vec sdf);
RcppExport SEXP _RiVIERA_lprPT_glass(SEXP pvalSEXP, SEXP priSEXP, SEXP annSEXP, SEXP alphaSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP groupMatSEXP, SEXP sdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< vec >::type pri(priSEXP);
    Rcpp::traits::input_parameter< mat >::type ann(annSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< mat >::type groupMat(groupMatSEXP);
    Rcpp::traits::input_parameter< vec >::type sdf(sdfSEXP);
    rcpp_result_gen = Rcpp::wrap(lprPT_glass(pval, pri, ann, alpha, ann_w, ann_w0, ann_w0_mu, groupMat, sdf));
    return rcpp_result_gen;
END_RCPP
}
// lprPT_glass_givenPPA
double lprPT_glass_givenPPA(vec pval, vec pri, vec ppa, double alpha, vec ann_w, double ann_w0, double ann_w0_mu, mat groupMat, vec sdf);
RcppExport SEXP _RiVIERA_lprPT_glass_givenPPA(SEXP pvalSEXP, SEXP priSEXP, SEXP ppaSEXP, SEXP alphaSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP groupMatSEXP, SEXP sdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< vec >::type pri(priSEXP);
    Rcpp::traits::input_parameter< vec >::type ppa(ppaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< mat >::type groupMat(groupMatSEXP);
    Rcpp::traits::input_parameter< vec >::type sdf(sdfSEXP);
    rcpp_result_gen = Rcpp::wrap(lprPT_glass_givenPPA(pval, pri, ppa, alpha, ann_w, ann_w0, ann_w0_mu, groupMat, sdf));
    return rcpp_result_gen;
END_RCPP
}
// printAnn
void printAnn(std::vector< std::string > x, uvec activeGroupIdx, vec groupScores);
RcppExport SEXP _RiVIERA_printAnn(SEXP xSEXP, SEXP activeGroupIdxSEXP, SEXP groupScoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< std::string > >::type x(xSEXP);
    Rcpp::traits::input_parameter< uvec >::type activeGroupIdx(activeGroupIdxSEXP);
    Rcpp::traits::input_parameter< vec >::type groupScores(groupScoresSEXP);
    printAnn(x, activeGroupIdx, groupScores);
    return R_NilValue;
END_RCPP
}
// riviera_glass
List riviera_glass(vec pval, const mat& ann, mat groupMat, uvec targetGroups, std::vector<std::string> targetGroupNames, double minPval, double alpha0, double pri0, double epsilon, double thres, int max_iter, double minAlpha, double maxAlpha, double groupScoreThres, int verbose);
RcppExport SEXP _RiVIERA_riviera_glass(SEXP pvalSEXP, SEXP annSEXP, SEXP groupMatSEXP, SEXP targetGroupsSEXP, SEXP targetGroupNamesSEXP, SEXP minPvalSEXP, SEXP alpha0SEXP, SEXP pri0SEXP, SEXP epsilonSEXP, SEXP thresSEXP, SEXP max_iterSEXP, SEXP minAlphaSEXP, SEXP maxAlphaSEXP, SEXP groupScoreThresSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< mat >::type groupMat(groupMatSEXP);
    Rcpp::traits::input_parameter< uvec >::type targetGroups(targetGroupsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type targetGroupNames(targetGroupNamesSEXP);
    Rcpp::traits::input_parameter< double >::type minPval(minPvalSEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< double >::type pri0(pri0SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type minAlpha(minAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type maxAlpha(maxAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type groupScoreThres(groupScoreThresSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(riviera_glass(pval, ann, groupMat, targetGroups, targetGroupNames, minPval, alpha0, pri0, epsilon, thres, max_iter, minAlpha, maxAlpha, groupScoreThres, verbose));
    return rcpp_result_gen;
END_RCPP
}
// calGradPT_ridge
vec calGradPT_ridge(const vec& pri, const vec& ppa, const mat& ann, const vec& ann_w, const vec& ann_w_mu, double ann_w0, double ann_w0_mu);
RcppExport SEXP _RiVIERA_calGradPT_ridge(SEXP priSEXP, SEXP ppaSEXP, SEXP annSEXP, SEXP ann_wSEXP, SEXP ann_w_muSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type pri(priSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ppa(ppaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    rcpp_result_gen = Rcpp::wrap(calGradPT_ridge(pri, ppa, ann, ann_w, ann_w_mu, ann_w0, ann_w0_mu));
    return rcpp_result_gen;
END_RCPP
}
// calHessPT_ridge
mat calHessPT_ridge(const vec& pri, const mat& ann);
RcppExport SEXP _RiVIERA_calHessPT_ridge(SEXP priSEXP, SEXP annSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type pri(priSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    rcpp_result_gen = Rcpp::wrap(calHessPT_ridge(pri, ann));
    return rcpp_result_gen;
END_RCPP
}
// lprPT_ridge
double lprPT_ridge(vec pval, vec pri, mat ann, double alpha, vec ann_w, vec ann_w_mu, double ann_w0, double ann_w0_mu);
RcppExport SEXP _RiVIERA_lprPT_ridge(SEXP pvalSEXP, SEXP priSEXP, SEXP annSEXP, SEXP alphaSEXP, SEXP ann_wSEXP, SEXP ann_w_muSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< vec >::type pri(priSEXP);
    Rcpp::traits::input_parameter< mat >::type ann(annSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    rcpp_result_gen = Rcpp::wrap(lprPT_ridge(pval, pri, ann, alpha, ann_w, ann_w_mu, ann_w0, ann_w0_mu));
    return rcpp_result_gen;
END_RCPP
}
// lprPT_givenPPA_ridge
double lprPT_givenPPA_ridge(vec pval, vec pri, vec ppa, mat ann, double alpha, vec ann_w, vec ann_w_mu, double ann_w0, double ann_w0_mu);
RcppExport SEXP _RiVIERA_lprPT_givenPPA_ridge(SEXP pvalSEXP, SEXP priSEXP, SEXP ppaSEXP, SEXP annSEXP, SEXP alphaSEXP, SEXP ann_wSEXP, SEXP ann_w_muSEXP, SEXP ann_w0SEXP, SEXP ann_w0_muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< vec >::type pri(priSEXP);
    Rcpp::traits::input_parameter< vec >::type ppa(ppaSEXP);
    Rcpp::traits::input_parameter< mat >::type ann(annSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    rcpp_result_gen = Rcpp::wrap(lprPT_givenPPA_ridge(pval, pri, ppa, ann, alpha, ann_w, ann_w_mu, ann_w0, ann_w0_mu));
    return rcpp_result_gen;
END_RCPP
}
// riviera_ridge
List riviera_ridge(vec pval, const mat& ann, vec ann_w_mu, double minPval, double alpha0, double pri0, double epsilon, double thres, int max_iter, double minAlpha, double maxAlpha, bool fitAlpha);
RcppExport SEXP _RiVIERA_riviera_ridge(SEXP pvalSEXP, SEXP annSEXP, SEXP ann_w_muSEXP, SEXP minPvalSEXP, SEXP alpha0SEXP, SEXP pri0SEXP, SEXP epsilonSEXP, SEXP thresSEXP, SEXP max_iterSEXP, SEXP minAlphaSEXP, SEXP maxAlphaSEXP, SEXP fitAlphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type pval(pvalSEXP);
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w_mu(ann_w_muSEXP);
    Rcpp::traits::input_parameter< double >::type minPval(minPvalSEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< double >::type pri0(pri0SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type minAlpha(minAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type maxAlpha(maxAlphaSEXP);
    Rcpp::traits::input_parameter< bool >::type fitAlpha(fitAlphaSEXP);
    rcpp_result_gen = Rcpp::wrap(riviera_ridge(pval, ann, ann_w_mu, minPval, alpha0, pri0, epsilon, thres, max_iter, minAlpha, maxAlpha, fitAlpha));
    return rcpp_result_gen;
END_RCPP
}
// eval_lld
double eval_lld(vec zscore, mat ldmat, double s2);
RcppExport SEXP _RiVIERA_eval_lld(SEXP zscoreSEXP, SEXP ldmatSEXP, SEXP s2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type zscore(zscoreSEXP);
    Rcpp::traits::input_parameter< mat >::type ldmat(ldmatSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    rcpp_result_gen = Rcpp::wrap(eval_lld(zscore, ldmat, s2));
    return rcpp_result_gen;
END_RCPP
}
// eval_pri
double eval_pri(const mat& ann, const uvec& causalidx, const vec& ann_w, double ann_w0, bool logd);
RcppExport SEXP _RiVIERA_eval_pri(SEXP annSEXP, SEXP causalidxSEXP, SEXP ann_wSEXP, SEXP ann_w0SEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type ann(annSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type causalidx(causalidxSEXP);
    Rcpp::traits::input_parameter< const vec& >::type ann_w(ann_wSEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_pri(ann, causalidx, ann_w, ann_w0, logd));
    return rcpp_result_gen;
END_RCPP
}
// find_uniqueRows
umat find_uniqueRows(umat x);
RcppExport SEXP _RiVIERA_find_uniqueRows(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< umat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(find_uniqueRows(x));
    return rcpp_result_gen;
END_RCPP
}
// sampleConfigs
umat sampleConfigs(const vec& proposal, int sampleSize);
RcppExport SEXP _RiVIERA_sampleConfigs(SEXP proposalSEXP, SEXP sampleSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< int >::type sampleSize(sampleSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleConfigs(proposal, sampleSize));
    return rcpp_result_gen;
END_RCPP
}
// sampleAnnLambda
mat sampleAnnLambda(double wishart_df, mat ann_w_sigma0, vec ann_w);
RcppExport SEXP _RiVIERA_sampleAnnLambda(SEXP wishart_dfSEXP, SEXP ann_w_sigma0SEXP, SEXP ann_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type wishart_df(wishart_dfSEXP);
    Rcpp::traits::input_parameter< mat >::type ann_w_sigma0(ann_w_sigma0SEXP);
    Rcpp::traits::input_parameter< vec >::type ann_w(ann_wSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleAnnLambda(wishart_df, ann_w_sigma0, ann_w));
    return rcpp_result_gen;
END_RCPP
}
// sampleAnnBiasTau
double sampleAnnBiasTau(double ann_w0, double ann_w0_mu, double alpha, double beta);
RcppExport SEXP _RiVIERA_sampleAnnBiasTau(SEXP ann_w0SEXP, SEXP ann_w0_muSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type ann_w0(ann_w0SEXP);
    Rcpp::traits::input_parameter< double >::type ann_w0_mu(ann_w0_muSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleAnnBiasTau(ann_w0, ann_w0_mu, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RiVIERA_dbetaPT", (DL_FUNC) &_RiVIERA_dbetaPT, 2},
    {"_RiVIERA_inferPri", (DL_FUNC) &_RiVIERA_inferPri, 3},
    {"_RiVIERA_inferPriByLoci", (DL_FUNC) &_RiVIERA_inferPriByLoci, 4},
    {"_RiVIERA_rwishart", (DL_FUNC) &_RiVIERA_rwishart, 2},
    {"_RiVIERA_isPosDef", (DL_FUNC) &_RiVIERA_isPosDef, 1},
    {"_RiVIERA_dmvnrm", (DL_FUNC) &_RiVIERA_dmvnrm, 4},
    {"_RiVIERA_dmvnrm_lambda", (DL_FUNC) &_RiVIERA_dmvnrm_lambda, 4},
    {"_RiVIERA_fastInv", (DL_FUNC) &_RiVIERA_fastInv, 1},
    {"_RiVIERA_duvnrm", (DL_FUNC) &_RiVIERA_duvnrm, 3},
    {"_RiVIERA_createLocusCursor", (DL_FUNC) &_RiVIERA_createLocusCursor, 2},
    {"_RiVIERA_normalizeZscoreByLoci", (DL_FUNC) &_RiVIERA_normalizeZscoreByLoci, 2},
    {"_RiVIERA_estimateParByEM", (DL_FUNC) &_RiVIERA_estimateParByEM, 6},
    {"_RiVIERA_calGradW_l", (DL_FUNC) &_RiVIERA_calGradW_l, 9},
    {"_RiVIERA_calGradW", (DL_FUNC) &_RiVIERA_calGradW, 10},
    {"_RiVIERA_calGrad_dEdV", (DL_FUNC) &_RiVIERA_calGrad_dEdV, 2},
    {"_RiVIERA_calGrad_dEdS", (DL_FUNC) &_RiVIERA_calGrad_dEdS, 3},
    {"_RiVIERA_calGradS_l", (DL_FUNC) &_RiVIERA_calGradS_l, 5},
    {"_RiVIERA_calGradS", (DL_FUNC) &_RiVIERA_calGradS, 6},
    {"_RiVIERA_hmcS", (DL_FUNC) &_RiVIERA_hmcS, 11},
    {"_RiVIERA_hmcW", (DL_FUNC) &_RiVIERA_hmcW, 13},
    {"_RiVIERA_inferConfigs_imp", (DL_FUNC) &_RiVIERA_inferConfigs_imp, 9},
    {"_RiVIERA_inferConfigs_nbh", (DL_FUNC) &_RiVIERA_inferConfigs_nbh, 8},
    {"_RiVIERA_inferConfigsAllLoci", (DL_FUNC) &_RiVIERA_inferConfigsAllLoci, 11},
    {"_RiVIERA_inferPPA", (DL_FUNC) &_RiVIERA_inferPPA, 6},
    {"_RiVIERA_inferPPAByLoci", (DL_FUNC) &_RiVIERA_inferPPAByLoci, 6},
    {"_RiVIERA_initParams", (DL_FUNC) &_RiVIERA_initParams, 7},
    {"_RiVIERA_lprW", (DL_FUNC) &_RiVIERA_lprW, 10},
    {"_RiVIERA_lprS_l", (DL_FUNC) &_RiVIERA_lprS_l, 5},
    {"_RiVIERA_nbh", (DL_FUNC) &_RiVIERA_nbh, 2},
    {"_RiVIERA_riviera_fmap", (DL_FUNC) &_RiVIERA_riviera_fmap, 22},
    {"_RiVIERA_calGradPT_glass", (DL_FUNC) &_RiVIERA_calGradPT_glass, 8},
    {"_RiVIERA_calHessPT_glass", (DL_FUNC) &_RiVIERA_calHessPT_glass, 5},
    {"_RiVIERA_lprPT_glass", (DL_FUNC) &_RiVIERA_lprPT_glass, 9},
    {"_RiVIERA_lprPT_glass_givenPPA", (DL_FUNC) &_RiVIERA_lprPT_glass_givenPPA, 9},
    {"_RiVIERA_printAnn", (DL_FUNC) &_RiVIERA_printAnn, 3},
    {"_RiVIERA_riviera_glass", (DL_FUNC) &_RiVIERA_riviera_glass, 15},
    {"_RiVIERA_calGradPT_ridge", (DL_FUNC) &_RiVIERA_calGradPT_ridge, 7},
    {"_RiVIERA_calHessPT_ridge", (DL_FUNC) &_RiVIERA_calHessPT_ridge, 2},
    {"_RiVIERA_lprPT_ridge", (DL_FUNC) &_RiVIERA_lprPT_ridge, 8},
    {"_RiVIERA_lprPT_givenPPA_ridge", (DL_FUNC) &_RiVIERA_lprPT_givenPPA_ridge, 9},
    {"_RiVIERA_riviera_ridge", (DL_FUNC) &_RiVIERA_riviera_ridge, 12},
    {"_RiVIERA_eval_lld", (DL_FUNC) &_RiVIERA_eval_lld, 3},
    {"_RiVIERA_eval_pri", (DL_FUNC) &_RiVIERA_eval_pri, 5},
    {"_RiVIERA_find_uniqueRows", (DL_FUNC) &_RiVIERA_find_uniqueRows, 1},
    {"_RiVIERA_sampleConfigs", (DL_FUNC) &_RiVIERA_sampleConfigs, 2},
    {"_RiVIERA_sampleAnnLambda", (DL_FUNC) &_RiVIERA_sampleAnnLambda, 3},
    {"_RiVIERA_sampleAnnBiasTau", (DL_FUNC) &_RiVIERA_sampleAnnBiasTau, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_RiVIERA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
