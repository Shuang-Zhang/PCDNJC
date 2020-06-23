#include "cdcparallel.h"
#include <chrono>
#include<algorithm>

void cdcparallel(Matrix &x,double &lam,Ithist &ithist,
                const Matrix& X,const Matrix &y,Opts opts){
    int n = X.Size(0);
    int p = X.Size(1);
    Matrix Xty(X);
    Xty.Transpose().Multiply(y);
    double linf = norm(Xty,NormType::inf);
    x.Zeros({p,1});
    // if(opts.var.Size() == 0){
    //     opts.N = 100;
    //     opts.Lmax = 1;
    //     opts.Lmin = 1e-10;
    //     opts.MaxIt = 5;
    //     opts.mu = 0.5;
    //     Matrix temp(y);
    //     //opts.del = norm(temp.Subtract(ye));
    // }
    double cnst = 0;;
    if(opts.method == "l0"){
        std::cout<<"Parallel l_0 model is running ..."<<std::endl;
        cnst = linf*linf /2;
        opts.gamma = 0;
    }else if(opts.method == "l1"){
        std::cout << "Parallel l_1 model is running ..."<<std::endl;
        cnst = linf;
        opts.gamma = 0;
    }else if(opts.method == "scad"){
        std::cout << "Parallel scad model is running ..."<<std::endl;
        cnst = linf;
        opts.gamma = 3.7;
    }else if(opts.method == "mcp"){
        std::cout <<"Parallel mcp model is running ..."<<std::endl;
        cnst = linf;
        opts.gamma = 1.5;
    }else{
        std::cout<<"Undefined Penaty"<<std::endl;
    }
    Matrix LamTemp;
    LamTemp.Linspace(log(opts.Lmax),log(opts.Lmin),opts.N).Exp();
    Matrix Lam;
    LamTemp.Get(Lam,1,LamTemp.Size(),0).Transpose().Multiply(cnst);
    ithist.Lam = Lam;
    Matrix xini = x;
    Matrix resini = y;
    int LamSize = Lam.Size();

    double res = 0;
    //ithist.res.Empty({LamSize,1});
    //ithist.bic.Empty({LamSize,1});
    //ithist.as.Empty({opts.numcore,LamSize});
    //ithist.it.Empty({opts.numcore,LamSize});
    //ithist.t.Empty({opts.numcore,LamSize});
    //ithist.x.Empty({x.Size(0),LamSize});
    for(int k=0;k<LamSize;++k){
        opts.lambda = Lam.Get(k);
        Matrix Tcore;
        Tcore.Zeros({opts.numcore,1});
        Matrix Kcore = Tcore;
        Matrix Itcore = Tcore;
        Matrix Procore;
        Procore.Zeros({n,opts.numcore});
        
#pragma omp parallel for num_threads(opts.thread)
        for(int nk = 0;nk<opts.numcore;++nk){
            auto start = std::chrono::high_resolution_clock::now();
            //由于c数组起始下标为0，Matlab数组起始下标为1,所以这里进行了调整
            std::vector<int> vecidxnk;
            Linspace(vecidxnk,nk*opts.pcore,(nk+1)*opts.pcore-1,opts.pcore);
            Matrix Xnk;
            X.Get(Xnk,vecidxnk,0);
            //这里多线程不能为opts.x0
            //xini.Get(opts.x0,vecidxnk,1);
            Matrix x0;
            xini.Get(x0,vecidxnk,1);
            Matrix ynk(Xnk);
            ynk.Multiply(x0).Add(resini);
            Matrix xnk,Xnkxnk;
            int snk,itnk;
            coordes(xnk,snk,itnk,Xnkxnk,Xnk,ynk,opts,x0);
            //x(idxnk,:) = xnk;
            for(int i=0;i<xnk.Size(0);++i){
                for(int j=0;j<xnk.Size(1);++j){
                    x[vecidxnk[i]][j] = xnk[i][j];
                }
            }
            Kcore.Get(nk) = snk;
            Itcore.Get(nk) = itnk;
            for(int i=0;i<Xnkxnk.Size(0);++i){
                Procore[i][nk] = Xnkxnk.Get(i);
            }
            auto end                           = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;        
            Tcore.Get(nk) = diff.count();
        }
        xini = x;
        Matrix Xxini;
        Procore.Sum(Xxini,1); 
        resini.Array(y);
        resini.Subtract(Xxini);
        res = norm(resini);
        Matrix s;
        Kcore.Sum(s);        
        // for(int i=0;i<x.Size(0);++i){
        //     ithist.x[i][k] = x.Get(i);

        // }

        // for(int i=0;i<opts.numcore;++i){
        //     ithist.as[i][k] = Kcore.Get(i);
        //     ithist.it[i][k] = Itcore.Get(i);
        //     ithist.t[i][k] = Tcore.Get(i);
        // }
        //ithist.res.Get(k) = res;
        //ithist.bic.Get(k) = 0.5*res*res + log(p)*s.Get(0)/n;
        ithist.x.Concatenate(x,1);
        ithist.as.Concatenate(Kcore,1);
        ithist.it.Concatenate(Itcore,1);
        ithist.res.Concatenate(Matrix(std::vector<std::vector<double>>({{res}})),0);  
        ithist.t.Concatenate(Tcore,1);
        double bic = 0.5*res*res + log(p)*s.Get(0)/n;
        ithist.bic.Concatenate(Matrix(std::vector<std::vector<double>>({{bic}})),0);
        
        if(res <= opts.del){
            std::cout<<"Discrepancy Principle is satisfied"<<std::endl;
            lam = Lam.Get(k);
            ithist.x.Get(x,{k},0);
            break;
        }
        if(s.Get(0) > opts.mu *n){
            std::cout<<"Too many nonzeros in the reconstruction ...."<<std::endl; 
            break;
        }
    }
    if(res > opts.del){
        std::cout<<"Discrepancy principle is not satisfied and BIC is used"<<std::endl;
    
    //ii = find(ithist.bic == min(ithist.bic));
    Matrix min;
    ithist.bic.Min(min);
    std::vector<int> ii;
    for(int i=0;i<ithist.bic.Size();++i){
       if(ithist.bic.Get(i) == min.Get(0)){
           ii.push_back(i);
       }
    }
    ithist.x.Get(x,{ii[ii.size()-1]},0);
    lam = Lam.Get(ii.size()-1);
    }
}

void coordes(Matrix &x,int &s, int &it,Matrix &XAxA,
             const Matrix & X,const Matrix &y,const Opts& opts,const Matrix & x0){
    int n = X.Size(0);
    int p = X.Size(1);
    //x = opts.x0;
    x = x0;
    
    double lambda = opts.lambda;
    double gamma = opts.gamma;
    Matrix r(y);
    Matrix Xx(X);
    r.Subtract(Xx.Multiply(x));  // r=y-X*x
    s=find(x).size();
    double rechange = opts.stop +1;
    it = 0;
    while(rechange > opts.stop && it <= opts.MaxIt){
        XAxA.Zeros({n,1});
        it++;
        Matrix xold(x);
        for(int j =0; j<p ;++j){
            Matrix Xj;
            X.Get(Xj,{j},0);
            double xoldj = xold.Get(j);
            Matrix temp(Xj);
            double zj = temp.Transpose().Multiply(r).Get(0) + xoldj;
            double xj =0;
            if(opts.method == "l0"){
                if(abs(zj) > sqrt(2*lambda)){
                    xj = zj;               
                }else{
                    xj = 0;
                }
            }else if(opts.method == "l1"){
                xj = sign(zj)*std::max((abs(zj)-lambda),(double)0);
            }else if(opts.method == "scad"){
                if(abs(zj) > gamma *lambda){
                    xj=zj;
                }else if(2*lambda <abs(zj) && abs(zj) <= gamma * lambda){
                    xj = sign(zj)*std::max((abs(zj)-gamma * lambda)/(gamma -1),(double)0)*(gamma -1)/(gamma-2);
                }else{
                    xj = sign(zj)*std::max((abs(zj)-lambda),0.);
                }
            }else if(opts.method =="mcp"){
                if(abs(zj) > gamma *lambda){
                    xj = zj;
                }else{
                    xj = sign(zj) * std::max((abs(zj)-lambda),0.)*gamma/(gamma-1);
                }
            }else{
                std::cout<<"Undefined Penaty"<<std::endl;
            }
            x.Get(j) = xj;
            if(xj != 0){
                Matrix temp(Xj);
                temp.Multiply(xj);
                XAxA.Add(temp);
            }
            temp.Array(Xj);
            temp.Multiply(xj-xoldj);
            r.Subtract(temp);
        }
        s = find(x).size();
        Matrix temp(x);
        rechange = norm(temp.Subtract(xold))/std::max(1.,norm(xold));
    }

}

bool gendata(Matrix &X,Matrix &y,Matrix &xe,Matrix &Ae,
             int n,int p,int K,double sigma,double ratio,double seednum,double rho){
    std::cout<<"Data is generating ..."<<std::endl;
    srand(seednum);
    xe.Zeros({p,1});
    std::vector<int> q;
    Matrix::RandVector(q,p);
    std::vector<int> A(q.begin(),q.begin()+K);
    if(ratio != 0){
        Matrix vx;
        vx.Rand({K,1}).Multiply(ratio);
        Matrix max,min;
        vx.Subtract(vx.Min(min).Get(0));
        vx.Divide(vx.Max(max).Multiply(ratio).Get(0));
        Matrix randn;
        randn.RandN({K,1},seednum);
        for(int i=0;i<K;++i){
            xe.Get(A[i]) = std::pow(10,vx.Get(i))*sign(rand()-RAND_MAX/2);  //xe(A) = 10.^vx.*sign(randn(K,1))
        }
    }else{
        for(int i=0;i<K;++i){
            xe.Get(A[i])=  sign(rand()-RAND_MAX/2);
        }
           
    }
    X.RandN({n,p});
    Matrix Sigma;
    if(rho != 0){
        Sigma.Zeros({p,p});   
        for(int k=0;k<p;++k){
            for(int l=0;l<p;++l){
                Sigma[k][l] = pow(rho,abs(k-l));
            }
        }
        Matrix chol;
        X.Multiply(Sigma.Chol(chol));
    }
    Matrix temp1,temp2;
    normalize(temp1,temp2,X);   
    X.Array(temp1);       //X = normalize(X)

    Matrix ye(X);
    ye.Multiply(xe);
    y.Array(ye);
    Matrix drand;
    drand.RandN({n,1}).Multiply(sigma);
    y.Add(drand);

    //A = find(xe)
    A =find(xe);
    
    std::cout<<"Generating is done"<<std::endl;
    
    return 0;

}

