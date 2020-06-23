#include"utils.h"
#include <time.h>
std::vector<int> find(const Matrix &x){
    std::vector<int> result;
    int size = x.Size();
    for(int i=0;i<size;++i){
        if(x.Get(i) != 0){
          result.push_back(i);
        }
    }
    return result;
}

Matrix sign(const Matrix & x){
    int size=x.Size();
    Matrix result;
    result.Empty({x.Size(0),x.Size(1)});
    for(int i=0;i<size;++i){
        if(x.Get(i) > 0){
            result.Get(i) = 1;
        }else if(x.Get(i) < 0){
            result.Get(i) = -1;
        }
    }
    return result;
}

double sign(const double &x){
    if( x > 0 ){
        return 1;
    }else if( x<0 ){
        return -1;
    }else{
        return 0;
    }
}

bool Normalize(Matrix & sX,Matrix& d,const Matrix &X){
    //按列求l2范数存入d，然后对X按列除以该范数存入sX
    int n = X.Size(0);
    int p = X.Size(1);
    d.Zeros({p,1});
    sX=X;
    for(int k=0;k<p;++k){
       double norm =0;
       for(int i=0;i<n;++i){
           norm += X[k][i]*X[k][i];
       }
       norm = sqrt(norm);
       d.Get(k) = 1./norm;
       for(int i=0;i<n;++i){
           sX[k][i] = sX[k][i] /norm;
       }
    }

    return 0;
}

double norm(const Matrix &x1,const Matrix& x2,NormType methods){
    double result=0;
    if(x1.Size() != x2.Size()){
        std::cout<<"norm的两个矩阵大小不一致"<<std::endl;
        return result;
    }
    switch(methods){
      case NormType::inf:{
          for(int i=0;i<x1.Size();++i){
              if(fabs(x1.Get(i) - x2.Get(i)) > result){
                  result = fabs(x1.Get(i) - x2.Get(i)); 
              }
          }
          break;
      }
      case NormType::l2:{
          for(int i=0;i<x1.Size();++i){
            double temp = x1.Get(i) - x2.Get(i);
            result += temp*temp;
          }
          result = sqrt(result);
          break;
      }
      default:{
          std::cout<<"norm 函数 methods 输入不对"<<std::endl;
          break;
      }
    }
    return result;
}

double norm(const Matrix &x,NormType methods){
    double result=0;
    switch(methods){
      case NormType::inf:{
          for(int i=0;i<x.Size();++i){
              if(fabs(x.Get(i)) > result){
                  result = fabs(x.Get(i)); 
              }
          }
          break;
      }
      case NormType::l2:{
          for(int i=0;i<x.Size();++i){
            result += x.Get(i)*x.Get(i);
          }
          result = sqrt(result);
          break;
      }
      default:{
          std::cout<<"norm 函数 methods 输入不对"<<std::endl;
          break;
      }
    }
    return result;
}

void normalize(Matrix & sX,Matrix &d,const Matrix & X){
  int n = X.Size(0);
  int p = X.Size(1);
  d.Zeros({p,1});
  sX.Array(X);
  for(int k=0;k<p;++k){
       double norm =0;
       for(int i = 0;i<n;++i){
         norm += X[i][k] *X[i][k];
       }
       norm = sqrt(norm);
       d[k][0] = 1./norm;
       for(int i=0;i<n;++i){
         sX[i][k] = sX[i][k]/norm;
       }
  }
}

char getTokenFromFile(std::ifstream &file, char delim, char *value) {
        char c;
        int count = 0;

        while (file.get(c)){
            if (c == delim) {
                break;
            } else if (c == '\0' || c == '\n') {
                break;
            } else {
                if (c == ' ') {
                    continue;
                }
                value[count++] = c;
            }            
        }
        //跳过连续的空格
        while(file.get(c)){
          if(c != ' '){
              file.unget();
              break;
          }
        }
        
        value[count] = 0;
        return c;
    }

void char2numeric(const char *value, double &lfvalue) {
    if (value) {
        // sscanf(value, "%lf", &lfvalue);
        lfvalue = atof(value);
    }
}

bool load(Matrix &x,std::string strfile,int irow,int icol){
       std::ifstream file;
        file.open(strfile, std::ios::in);
        if (file.fail()) {
            std::cout << "Read():Error:Unable to open file %s" << strfile << std::endl;
            return 1;
        }
        std::vector<std::vector<double>> vecx;
        vecx.resize(irow, std::vector<double>(icol));
        char c;
        char value[32];
        double temp = 0;  
        int row   = 0;
        int col   = 0;
        while (file.get(c)) {
            file.unget();
            c = getTokenFromFile(file, ' ', value);
            if(value[0] != '\0'){
                char2numeric(value, temp);
                vecx[row][col] = temp;
                if(col+1 < icol){
                    col++;
                }else{
                    row++;
                    col=0;
                }
            }  
        }
        x.Array(vecx);
        return 0;


}

void parse_command_line(int argc, char **argv, DataParams &dataparams,Opts& opts)
{
    //默认参数
    // dataparams.dataset = 0; //0表示用Matlab生成数据集，否则这里重新生成
    dataparams.p = 4096;
    dataparams.n=1024;
    dataparams.sigma = 0.05;
    dataparams.K=32;
    dataparams.ratio = 1;
    dataparams.rho = 0;
    
    dataparams.seednum = 0;
    
    opts.N     = 50;          // you can set to 100 to improve accuracy         
    opts.Lmax  = 1;
    opts.Lmin  = 1e-5;
    opts.MaxIt = 20;   
    opts.stop = 1e-3;
    opts.method = "mcp";
    opts.respath = "../draw/";
    opts.numcore = 8;
    opts.thread = 8;

    // parse options
	for(int i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
            // case 'd':
            //     dataparams.dataset = atoi(argv[i]);
            //     break;
            case 'p':
                dataparams.p=atoi(argv[i]);
                break;
            case 'n':
                dataparams.n=atoi(argv[i]);
                break;
            case 'a':
                dataparams.sigma=atof(argv[i]);
                break;
            case 'K':
                dataparams.K=atoi(argv[i]);
                break;
            case 'R':
                dataparams.ratio=atof(argv[i]);
                break;
            case 'r':
                dataparams.rho=atof(argv[i]);
                break;
            case 's':
                dataparams.seednum=atof(argv[i]);
                break;
            case 'o':
                opts.method=std::string(argv[i]);
                break;
            case 'N':
                opts.N=atoi(argv[i]);
                break;
            case 'M':
                opts.Lmax=atof(argv[i]);
                break;
            case 'm':
                opts.Lmin=atof(argv[i]);
                break;
            case 'i':
                opts.MaxIt=atoi(argv[i]);
                break;
            case 't':
                opts.stop=atoi(argv[i]);
                break;
            case 'T':
                opts.thread = atoi(argv[i]);
                break;
            case 'c':
                opts.numcore = atoi(argv[i]);
                break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
                exit_with_help();
		}
        opts.mu    = 1/log(dataparams.p); 
        opts.del = 1*dataparams.sigma*sqrt(dataparams.n);  // Only  using  BIC by setting  opts.del = 0;
        //opts.ncore = dataparams.n;
	}
}

void exit_with_help(){
    std::cout<<
    "-p:length of solution"<<std::endl<<
    "-n:number of samples"<<std::endl<<
    "-K:number of nonzeros element"<<std::endl<<
    "-a:sigma ,noise voise variance"<<std::endl<<
    "-R:ratio,Dinamic range (10^ratio)"<<std::endl<<
    "-o:method,you can choose l0,l1,mcp,scad or all"<<std::endl<<
    "-r:rho,coorelation"<<std::endl<<
    "-s:seenum,seed number"<<std::endl<<
    "-N:opts.N you can set to 100 to improve accuracy"<<std::endl<<
    "-M:opts.Lmax"<<std::endl<<
    "-m:opts.Lmin"<<std::endl<<
    "-i:opts.maxIt"<<std::endl<<
    "-t:otps.stop"<<std::endl<<
    "-T:opts.thread,程序中使用的线程数"<<std::endl<<
    "-c:opts.numcore,注意opts.numcore值应大于等于opts.thread"<<std::endl;
    exit(1);
}

void Linspace(std::vector<int> &result,double start, double stop, int num, bool endpoint){
        result.resize(num);
        double step = 0;
        if (endpoint == 1) {
            step = (stop - start) / (num - 1);
        } else {
            step = (stop - start) / num;
        }
        for (int i = 0; i < num; ++i) {
            result[i] = start + i * step;
        }

}