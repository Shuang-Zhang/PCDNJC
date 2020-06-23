#ifndef UTILS
#define UTILS

#include<string>
#include "matrix.h"

struct Opts{
  public:
    int numcore;
    double pcore;
    //int ncore;
    std::string method;   // mcp,l0,l1,mcp,scad
    std::string respath; 
    int N;                // you can set to 100 to improve accuracy 
    double Lmax;
    double Lmin;
    int MaxIt;
    double mu;
    double stop;
    double del;           // Only  using  BIC by setting  opts.del = 0;

    Matrix var;

    double lambda;
    double gamma;
    int thread;
};

struct DataParams{
    int dataset;

    int p;
    int n;
    double sigma;
    int K;
    double ratio;
    double rho;
    double seednum;
};

struct Ithist{
  Matrix Lam;
  Matrix x;
  Matrix as;
  Matrix it;

  Matrix res;
  Matrix bic;
  Matrix t;
};

enum class NormType{l2,inf};

double norm(const Matrix &x1,const Matrix& x2,NormType methods=NormType::l2,int axis=-1);

double norm(const Matrix &x,NormType methods = NormType::l2);

void normalize(Matrix & sX,Matrix &d,const Matrix & X);

std::vector<int> find(const Matrix &x);

Matrix sign(const Matrix & x);

double sign(const double &x);

bool load(Matrix &x,std::string strfile,int irow,int icol);

void parse_command_line(int argc, char **argv, DataParams &dataparams,Opts& opts);
void exit_with_help();

void Linspace(std::vector<int>& result,double start, double stop, int num, bool endpoint=1);
#endif