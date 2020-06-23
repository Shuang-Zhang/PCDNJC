// [[Rcpp::depends(Rcpp, RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <cmath>
#include<iostream>
#include<random>
#include <vector>       // std::vector
#include <iterator>
#include <string>

using namespace std;


double sign(double &x){
  if(x>0){
    return 1.0;
  }else if(x<0){
    return -1.0;
  }else{
    return 0;
  }
}


std::vector<Eigen::MatrixXd> normalize(Eigen::MatrixXd &x, Eigen::VectorXd &y, std::string nnn="sqrtn"){
  int m = x.cols();
  int n = x.rows();
  Eigen::VectorXd one = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd meanx = one.transpose() * x * (1/(double) n);
  for(int i = 0; i < n; ++i){
    x.row(i) = x.row(i) - meanx.transpose();
  }
  Eigen::VectorXd normx_inv(m);
  if(nnn=="sqrtn"){
    normx_inv = (one.transpose()*(x.unaryExpr([](double elem){return elem*elem;}))).unaryExpr([&n](double dummy){return sqrt((double) n)/sqrt(dummy);});
  }else{
    normx_inv = (one.transpose()*(x.unaryExpr([](double elem){return elem*elem;}))).unaryExpr([](double dummy){return 1.0/sqrt(dummy);});
  }
  double meany = y.sum()/(double) n;
  Eigen::VectorXd	yc = y.unaryExpr([&meany](double elem){return elem - meany;});
  //huge diag matrix for large p
  return {x*(normx_inv.asDiagonal()), yc, normx_inv};
}



//[[Rcpp::export]]
std::vector<Eigen::MatrixXd> gdata(const int &n=200, const int &p=1000, const int &K=10, double sigma=0.1, const double &ratio=1.0, const int &kind=1, double rho=1.0e-10, const int &seed=1, const bool &isnorm=true, std::string nnn="sqrtn", double snr=1.0e-10){
  //snr, rho 0.00000001 to capture no input
  //========================================================================
  // INPUTS:
  //     n   ---- number of samples
  //     p   ---- signal length
  //     K   ---- number of nonzero elements in the signal
  //  sigma  ---- noise variance
  //  ratio  ---- range of value in the signal (= 10^ratio)
  //   kind  ---- type of sample matirx
  //          1 for Random Gaussian with auto coorelation (small p)
  //          2 for correlated Random Gaussian  (large p)
  //    rho  ---- corr. coeff.
  // 	seed   ---- seed number
  // 	isnorm  ---- normalization or not
  //    nnn  ---- normalization index: \|x\|_2 = sqrt{n}  or  \|x\|_2 = 1, value in {"sqrtn", "one"}
  //    snr  ---- signal to noise ratio
  // OUTPUTS:
  //    btrue   ---- true signal
  //     A   ---- the support of btrue
  //     x,y,ye   ---- raw data
  //     nx,ny  ---- normalized data
  //     d   ---- vector makes raw data normalized
  //------------------------------------------------------------------------
  
  static std::default_random_engine e{seed};
  static normal_distribution<double> normdis(0,1);
  Eigen::VectorXd btrue = Eigen::VectorXd::Zero(p);
  std::vector<int> q;
  for(int i = 0; i < p; ++i){
    q.push_back(i);
  }
  std::random_shuffle(q.begin(),q.end());
  //Not decent way to get first K elem.
  std::vector<int> A;
  for(int i=0; i < K; ++i){
    A.push_back(q[i]);
  }
  if(ratio == 0){
    for(int j = 0; j < K; ++j){
      btrue(A[j]) = (normdis(e)>0)?1.0:-1.0;
    }
  }else{
    Eigen::VectorXd vx = Eigen::VectorXd::Zero(n).unaryExpr([&ratio](double dummy){return ratio*normdis(e);});
    double vx_min = vx.minCoeff();
    Eigen::VectorXd vx_tmp = vx.unaryExpr([&vx_min](double dummy){return dummy-vx_min;});
    double vx_max = vx_tmp.maxCoeff();
    vx = vx_tmp.unaryExpr([&ratio,vx_max](double dummy){return dummy*ratio/vx_max;});
    for(int j = 0; j < K; ++j){
      btrue(A[j]) = pow(10, vx(j))*((normdis(e)>0)?1.0:-1.0);
    }
  }
  Eigen::MatrixXd x(n,p);
  if(kind == 1){
    x = Eigen::MatrixXd::Zero(n,p).unaryExpr([](double dummy){return normdis(e);});
    if(rho == 1.0e-10){
      rho = 0;
    }
    if(rho!=0){
      Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(p,p);
      for(int i = 0; i < p; ++i){
        for(int j = 0; j < p; ++j){
          Sigma(i, j) = pow(rho, fabs(i-j));
        }
      }
      Eigen::LLT<Eigen::MatrixXd> lltOfSigma(Sigma);
      x = x * lltOfSigma.matrixL();
      //Huge matrix should be released.
    }
  }else{
    if(rho==1.0e-10){rho=0.2;}
    x = Eigen::MatrixXd::Zero(n,p).unaryExpr([](double dummy){return normdis(e);});
    Eigen::VectorXd tmp(n);
    for(int j = 1; j < (p-1); ++j){
      x.col(j) = x.col(j) + rho * x.col(j+1) + rho * x.col(j-1);
    }
  }
  Eigen::VectorXd ye = x * btrue;
  Eigen::VectorXd noise = Eigen::VectorXd::Zero(n).unaryExpr([](double dummy){return normdis(e);});
  Eigen::VectorXd y(n);
  if(snr==1.0e-10){
    y = ye + sigma * noise;
  }else{
    sigma = (ye.array() - ye.array().mean()).square().sum()/((n-1)*sqrt(snr));
    y = ye + sigma * noise;
  }
  Eigen::MatrixXd nx(n,p);
  Eigen::VectorXd ny(n);
  Eigen::VectorXd d(n);
  if(isnorm){
    std::vector<Eigen::MatrixXd> ndat;
    if(nnn == "sqrtn"){
      ndat = normalize(x, y, "sqrtn");
    }else{
      ndat = normalize(x, y, "one");
    }
    //y is not empty so ny exist.
    nx = ndat[0];
    ny = ndat[1];
    d = ndat[2];
  }else{
    nx = x;
    ny = y;
  }
  Eigen::VectorXd A_eigenvec(K);
  for(int i = 0; i < K; ++i){
    //C++ index start from 0
    A_eigenvec(i) = A[i]+1;
  }
  return {btrue, A_eigenvec, x, y, ye, nx, ny, d};
}


Eigen::VectorXd soft(Eigen::VectorXd &z, double lambda=1){
  return z.unaryExpr([&lambda](double elem){if(fabs(elem)>lambda){return sign(elem)*(fabs(elem)-lambda);}else{return 0.0;}});
}


Eigen::VectorXd thresh(Eigen::VectorXd &z, double lambda=1, double gamma=3, std::string pen="mcp"){
  //pen in "mcp","scad","lasso"
  if(lambda<0){Rcpp::Rcout << "lambda must be '>= 0'"<< std::endl;}
  if((gamma<=1) && (pen=="mcp")){Rcpp::Rcout <<"gamma for mcp must be '> 1'; gamma=3 is recommended"<< std::endl;}
  if((gamma<=2) && (pen=="scad")){Rcpp::Rcout <<"gamma for scad must be '> 2'; gamma=3.7 is recommended"<< std::endl;}
  Eigen::VectorXd th = Eigen::VectorXd::Zero(z.size());
  if(pen=="lasso"){
    th = soft(z, lambda);
  }else if(pen=="mcp"){
    Eigen::VectorXd zz_tmp(z.size());
    Eigen::VectorXi idx_tmp(z.size());
    int idx_index = 0;
    for(int i = 0; i < z.size(); ++i){
      if(fabs(z(i))<=(gamma*lambda)){
        idx_tmp(idx_index) = i;
        zz_tmp(idx_index) = z(i);
        ++idx_index;
      }else{th(i) = z(i);}
    }
    Eigen::VectorXd zz = zz_tmp.head(idx_index);
    Eigen::VectorXi idx = idx_tmp.head(idx_index);
    Eigen::VectorXd th_tmp = soft(zz, lambda)*(1/(1-1/gamma));
    for(int i = 0; i < idx.size(); ++i){
      th(idx(i)) = th_tmp(i);
    }
  }else if(pen=="scad"){
    Eigen::VectorXd zz_tmp1(z.size());
    Eigen::VectorXd zz_tmp2(z.size());
    Eigen::VectorXi idx_tmp1(z.size());
    Eigen::VectorXi idx_tmp2(z.size());
    int idx_index1 = 0;
    int idx_index2 = 0;
    for(int i = 0; i < z.size(); ++i){
      if(fabs(z(i))<=(2*lambda)){
        idx_tmp1(idx_index1) = i;
        zz_tmp1(idx_index1) = z(i);
        ++idx_index1;
      }else if(fabs(z(i))>(gamma*lambda)){
        th(i) = z(i);
      }else{
        idx_tmp2(idx_index2) = i;
        zz_tmp2(idx_index2) = z(i);
        ++idx_index2;
      }
    }
    Eigen::VectorXd zz1 = zz_tmp1.head(idx_index1);
    Eigen::VectorXd zz2 = zz_tmp2.head(idx_index2);
    Eigen::VectorXi idx1 = idx_tmp1.head(idx_index1);
    Eigen::VectorXi idx2 = idx_tmp2.head(idx_index2);
    Eigen::VectorXd th_tmp1 = soft(zz1, lambda);
    for(int i = 0; i < idx1.size(); ++i){
      th(idx1(i)) = th_tmp1(i);
    }
    Eigen::VectorXd th_tmp2 = soft(zz2, gamma*lambda/(gamma-1))*(1/(1-1/(gamma-1)));
    for(int i = 0; i < idx2.size(); ++i){
      th(idx2(i)) = th_tmp2(i);
    }
  }else{Rcpp::Rcout << "Expected penalty type!!!"<< std::endl;}
  return th;
}


std::vector<Eigen::VectorXd> cdj(Eigen::MatrixXd &x, Eigen::VectorXd &y, Eigen::VectorXd &beta0, std::string pen="mcp", std::string nnn="sqrtn", double la=1.0, double ga=3.0, double mu=0.5, int imax=5, double tol=1.0e-5){
  // for later calculate.
  double n = (double) x.rows();
  Eigen::VectorXd ebeta = beta0;
  Eigen::VectorXd tbeta = beta0;
  double as = ebeta.unaryExpr([](double elem){if(elem != 0){return 1.0;}else{return 0.0;}}).sum();
  double it = 0;
  double mun = mu * n;
  double change = 0.0;
  Eigen::VectorXd vec_as(1);
  Eigen::VectorXd vec_it(1);
  while((as<mun)&&(it<imax)){
    if(nnn=="sqrtn"){
      tbeta = beta0 + (1/n)*x.transpose()*(y-x*beta0);
    }else{
      tbeta = beta0 + x.transpose()*(y-x*beta0);
    }
    ebeta = thresh(tbeta,la,ga,pen);
    as = ebeta.unaryExpr([](double elem){if(elem != 0){return 1.0;}else{return 0.0;}}).sum();
    change = sqrt((beta0-ebeta).unaryExpr([](double elem){return elem*elem;}).sum());
    if(change<tol){
      break;
    }
    beta0 = ebeta;
    ++it;
  }
  vec_as(0) = as;
  vec_it(0) = it;
  return {ebeta, vec_as, vec_it};
}


std::vector<Eigen::VectorXd> lamfun(Eigen::MatrixXd &x, Eigen::VectorXd &y, int N=100, std::string nnn="sqrtn", double Lmin=1.0e-10){
  //pen in "mcp","scad"
  int n = x.rows();
  int p = x.cols();
  Eigen::VectorXd linf_vec(1);
  Eigen::VectorXd cnst_vec(1);
  if(nnn == "sqrtn"){
    linf_vec(0) = ((x.transpose()*y*(1/(double) n)).unaryExpr([](double elem){return fabs(elem);})).maxCoeff();
  }else{
    linf_vec(0) = ((x.transpose()*y).unaryExpr([](double elem){return fabs(elem);})).maxCoeff();
  }
  if(Lmin==1.0e-10){
    if(n>p){
      Lmin = 1.0e-3;
    }else{
      Lmin = 1.0e-5; 
    }
  }
  Eigen::VectorXd Lam(N);
  double step = (log(Lmin)-log(1))/(N-1);
  Lam(0) = log(1);
  for(int i = 1; i < N; ++i){
    Lam(i) = log(1)+i*step;
  }
  for(int i = 0; i < N; ++i){
    Lam(i) = exp(Lam(i));
  }
  cnst_vec = linf_vec;
  Lam = Lam*cnst_vec(0);
  return {linf_vec, cnst_vec, Lam};
}

//[[Rcpp::export]]
std::vector<Eigen::MatrixXd> cdnjc(Eigen::MatrixXd &x, Eigen::VectorXd &y, std::string pen="mcp", std::string nnn="sqrtn", int N=100, double Lmin=1.0e-5, double ga=1.0e-10, double mu=0.5, int imax=1, double tol=1.0e-5){
  //1.0e-10 capture for null.
  if(ga==1.0e-10){
    if(pen=="mcp"){ga=3;}
    else if(pen=="scad"){ga=3.7;}
  }
  //for later calculation
  double n = (double) x.rows();
  int p = x.cols();
  std::vector<Eigen::VectorXd> Laa = lamfun(x,y,N,nnn,Lmin);
  Eigen::VectorXd La = Laa[2];
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd As = Eigen::VectorXd::Zero(N);
  Eigen::MatrixXd Beta = Eigen::MatrixXd::Zero(p, N);
  Eigen::VectorXd Bic = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd It = Eigen::VectorXd::Zero(N);
  double la;
  std::vector<Eigen::VectorXd> tmp;
  Eigen::VectorXd bela(p);
  double asla, itla, res;
  for(int k = 0; k < N; ++k){
    la = La(k);
    tmp = cdj(x,y,beta0,pen,nnn,la,ga,mu,imax,tol);
    bela = tmp[0];
    asla = tmp[1](0);
    itla = tmp[2](0);
    Beta.col(k) = bela;
    As(k) = asla;
    It(k) = itla;
    res = sqrt((x*bela-y).unaryExpr([](double elem){return elem*elem;}).sum());
    if(nnn=="sqrtn"){
      Bic(k) = n*log(res*res/n) + log(n)*log(p)*asla;
    }else{
      Bic(k) = res*res + log(log(n))*log(p)*asla;
    }
    beta0 = bela;
  }
  int id = 0;
  double bic_min = Bic.minCoeff();
  for(int i = 0; i < N; ++i){
    if(Bic(i) == bic_min){
      id = i;
    }
  }
  Eigen::VectorXd lahat_vec(1);
  Eigen::VectorXd id_vec(1);
  lahat_vec(0) = La(id);
  id_vec(0) = id;
  return {Beta.col(id), lahat_vec, id_vec, As, Beta, Bic, It, La};
}







