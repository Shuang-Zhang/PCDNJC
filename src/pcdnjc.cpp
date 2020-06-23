#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "matrix.h"
#include "cdcparallel.h"
#include <chrono>


#include "openblas/cblas.h"
#include "openblas/lapacke.h"

#include <fstream>
#include <ostream>
#include <cstdio>



using namespace std::chrono;
using namespace std;


void test_parallel(const Matrix &xe,const Matrix&X,const Matrix&y,Opts opts,
                   const Matrix &Ae,const DataParams & dataParams);

// [[Rcpp::export]]
int pcdnjc(int p=4096, int n=1024, int K=32, double sigma=0.05, double ratio=1.0, int seednum=0, double rho=0.0, 
        double mu=0.1202, double del=1.6, int thread=8, int numcore=8, const char* method="all", int N=50, 
        double Lmax=1.0, double Lmin=1e-5, int MaxIt=20, double stop = 1e-3, const char* respath="draw/"){
    //设置openblas线程数为1
    openblas_set_num_threads(1);
    Matrix X,y,xe,Ae;
    Opts opts;
    DataParams pp;
    //parse_command_line(argc,argv, p,opts);
    pp.p = p;
    //pp.dataset = dataset;
    pp.n=n;
    pp.K=K;
    pp.sigma=sigma;
    pp.ratio=ratio;
    pp.seednum=seednum;
    pp.rho=rho;
    opts.mu    = mu; 
    opts.del = del;
    opts.thread=thread;
    opts.numcore=numcore;
    opts.method=method;
    opts.respath=respath;
    opts.N = N;
    opts.Lmax = Lmax;
    opts.Lmin = Lmin;
    opts.MaxIt = MaxIt;
    opts.stop = stop;
    // if(pp.dataset == 0){
    //     cout<<"默认使用旧有的由Matlab固定seednum参数为0生成的数据集"<<endl;
    //     cout<<"若要重新生成数据集，请使用 -d 参数并设置非零值"<<endl;
    //     load(X,"data/X.txt",1024,4096);
    //     load(y,"data/y.txt",1024,1);
    //     load(xe,"data/xe.txt",4096,1);
    //     load(Ae,"data/Ae.txt",32,1);
    //     pp.p = 4096;
    //     pp.n=1024;
    //     opts.mu    = 1/log(pp.p); 
    //     opts.del = 1*pp.sigma*sqrt(pp.n);  // Only  using  BIC by setting  opts.del = 0;
    // }else{
        gendata(X,y,xe,Ae,pp.n,pp.p,pp.K,pp.sigma,pp.ratio,pp.seednum,pp.rho);
    // }
    if(opts.method=="all"){
        std::vector<std::string> strmethod={"mcp","l0","l1","scad"};
        for(auto &iter :strmethod){
            opts.method = iter;
            test_parallel(xe,X,y,opts,Ae,pp);
        }
    }else{
        test_parallel(xe,X,y,opts,Ae,pp);
    }
    return 0;
    
}

void test_parallel(const Matrix &xe,const Matrix&X,const Matrix&y,Opts opts,
                   const Matrix & Ae,const DataParams &dataParams){
    std::cout<<"--------"<<opts.method<<"---------"<<std::endl;
    int maxthread = opts.thread;
    
    std::string strfile = opts.respath + opts.method+to_string(dataParams.dataset)+".txt";
    remove(strfile.c_str());
    //记录时间等结果，将用于在Matlab中绘图
    std::ofstream outfile;
    outfile.open(strfile);
    outfile <<"p:"<<endl<<dataParams.p<<endl;
    outfile <<"K:"<<endl<<dataParams.K<<endl;
    outfile << "xe:"<<endl<<xe<<endl;
    outfile << "Ae:"<<endl<<Ae<<endl;
    for(int i=1;i<=maxthread;i*=2){
        outfile << "thread:" << endl<< i << endl;
        opts.pcore = dataParams.p/opts.numcore;
        opts.thread = i;
        auto start = std::chrono::high_resolution_clock::now();
        std::cout<<"opts.numcore:"<<opts.numcore<<::endl;
        std::cout<<"线程数为:"<<opts.thread<<std::endl;
        Matrix x;
        double lam;
        Ithist ithist;
        cdcparallel(x,lam,ithist,X,y,opts);
        auto end                           = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        std::cout<<"总运行时间："<<diff.count()<<std::endl;
        Matrix error(x);
        error.Subtract(xe);
        double err_l2 = norm(error)/norm(xe);
        double err_linf = norm(error,NormType::inf);
        std::cout<<"err_l2:"<<err_l2<<std::endl;
        std::cout<<"err_linf:"<<err_linf<<std::endl<<endl;        
        
        if(i == 1){
            outfile << "x:"<<endl<<x<<endl;
            outfile << "ithist.res:"<<endl<<ithist.res<<endl;
            outfile <<"ithist.x:"<<endl<<ithist.x<<endl;
            outfile << "ithist.as:" <<endl<<ithist.as<<endl;
            outfile << "lam:"<<endl<< lam<<endl;
            outfile << "err_l2:"<<endl<<err_l2<<endl;
            outfile << "err_linf:"<<endl<<err_linf<<endl;
            
        }
        outfile <<"total_time:"<<endl<<diff.count()<<endl;  
    }
    outfile.close();
}
