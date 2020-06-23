// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// gdata
std::vector<Eigen::MatrixXd> gdata(const int& n, const int& p, const int& K, double sigma, const double& ratio, const int& kind, double rho, const int& seed, const bool& isnorm, std::string nnn, double snr);
RcppExport SEXP _CDJNC_gdata(SEXP nSEXP, SEXP pSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP ratioSEXP, SEXP kindSEXP, SEXP rhoSEXP, SEXP seedSEXP, SEXP isnormSEXP, SEXP nnnSEXP, SEXP snrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type ratio(ratioSEXP);
    Rcpp::traits::input_parameter< const int& >::type kind(kindSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type isnorm(isnormSEXP);
    Rcpp::traits::input_parameter< std::string >::type nnn(nnnSEXP);
    Rcpp::traits::input_parameter< double >::type snr(snrSEXP);
    rcpp_result_gen = Rcpp::wrap(gdata(n, p, K, sigma, ratio, kind, rho, seed, isnorm, nnn, snr));
    return rcpp_result_gen;
END_RCPP
}
// cdnjc
std::vector<Eigen::MatrixXd> cdnjc(Eigen::MatrixXd& x, Eigen::VectorXd& y, std::string pen, std::string nnn, int N, double Lmin, double ga, double mu, int imax, double tol);
RcppExport SEXP _CDJNC_cdnjc(SEXP xSEXP, SEXP ySEXP, SEXP penSEXP, SEXP nnnSEXP, SEXP NSEXP, SEXP LminSEXP, SEXP gaSEXP, SEXP muSEXP, SEXP imaxSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::string >::type pen(penSEXP);
    Rcpp::traits::input_parameter< std::string >::type nnn(nnnSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< double >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type imax(imaxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cdnjc(x, y, pen, nnn, N, Lmin, ga, mu, imax, tol));
    return rcpp_result_gen;
END_RCPP
}
// pcdnjc
int pcdnjc(int p, int dataset, int n, int K, double sigma, double ratio, int seednum, double rho, double mu, double del, int thread, int numcore, const char* method, int N, double Lmax, double Lmin, int MaxIt, double stop, const char* respath);
RcppExport SEXP _CDJNC_pcdnjc(SEXP pSEXP, SEXP datasetSEXP, SEXP nSEXP, SEXP KSEXP, SEXP sigmaSEXP, SEXP ratioSEXP, SEXP seednumSEXP, SEXP rhoSEXP, SEXP muSEXP, SEXP delSEXP, SEXP threadSEXP, SEXP numcoreSEXP, SEXP methodSEXP, SEXP NSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP MaxItSEXP, SEXP stopSEXP, SEXP respathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    Rcpp::traits::input_parameter< int >::type seednum(seednumSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type del(delSEXP);
    Rcpp::traits::input_parameter< int >::type thread(threadSEXP);
    Rcpp::traits::input_parameter< int >::type numcore(numcoreSEXP);
    Rcpp::traits::input_parameter< const char* >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< double >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type MaxIt(MaxItSEXP);
    Rcpp::traits::input_parameter< double >::type stop(stopSEXP);
    Rcpp::traits::input_parameter< const char* >::type respath(respathSEXP);
    rcpp_result_gen = Rcpp::wrap(pcdnjc(p, dataset, n, K, sigma, ratio, seednum, rho, mu, del, thread, numcore, method, N, Lmax, Lmin, MaxIt, stop, respath));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CDJNC_gdata", (DL_FUNC) &_CDJNC_gdata, 11},
    {"_CDJNC_cdnjc", (DL_FUNC) &_CDJNC_cdnjc, 10},
    {"_CDJNC_pcdnjc", (DL_FUNC) &_CDJNC_pcdnjc, 19},
    {NULL, NULL, 0}
};

RcppExport void R_init_CDJNC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
