#ifndef CDPARALLEL
#define CDPARALLEL

#include"matrix.h"
#include"utils.h"


bool gendata(Matrix &X,Matrix &y,Matrix &xe,Matrix &Ae,
             int n,int p,int K,double sigma,double ratio,double seednum,double rho);

void cdcparallel(Matrix &x,double &lam,Ithist &ithist,
                const Matrix& X,const Matrix &y,Opts opts);

void coordes(Matrix &x,int &s, int &it,Matrix &XAxA,
             const Matrix & X,const Matrix &y,const Opts& opts,const Matrix &x0);

#endif