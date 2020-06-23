# PCDNJC
R package "PCDNJC" for the Jacobi type CD algorithm with warmstarts in high-dimension sparse estimation, which is proposed by Jiao et al. (2020+).
# Installation
	# Requirs openblas envirnment for parallel computation.
	# Pure stand-alone version(Shuang-Zhang/CDNJC) can be used without openblas envirnment.
    # install.packages("devtools")
    library(devtools)
    install_github("Shuang-Zhang/PCDNJC")
   
# Example

    library(PCDNJC)
	# Stand-alone Version
	
	# Settings
	n=200; p=1000; K=10; sigma=0.1; ratio=1; kind=1; rho=0.1; seed=1; isnorm=TRUE; pen="scad"; nnn="sqrtn"; N=100; Lmin=1e-5; ga=4; mu=0.5; imax=1; tol=1e-5
	
	# Generate data
	data = gdata(n, p, K, sigma, ratio, kind, rho, seed, isnorm, nnn="sqrtn", snr=1.0e-10)
	
	# CDJ algorithm
	res_c = cdnjc(data[[6]],data[[7]],pen,nnn,N,Lmin,ga,mu,imax,tol) 
	
	# Parallel Version
	
	# Settings
	dataset=1 # 0 means using the prepared data, otherwise generate new data.
	method="all" # Choose from l0,l1,mcp,scad or all
	p=4096; n=1024; K=32; sigma=0.05; ratio=1.0; rho=0.0; seednum=0; Lmax=1.0; Lmin=1e-5; MaxIt=20; stop=1e-3; thread=8; numcore=8
	N=50 # Set to 100 to improve accuracy
	
	# PCDJ algorithm
	pcdnjc(p=p, dataset=dataset, n=n, K=K, sigma=sigma, ratio=ratio, seednum=seednum, rho=rho, thread=thread, numcore=numcore, method=method, N=N, Lmax=Lmax, Lmin=Lmin, MaxIt=MaxIt) #Result saved in draw filefolder
    
# References
Jiao Y. Coordinate Descent, Jacobi or Gauss-Siedel. Manuscript.

# Development
This R package is developed by Shuang Zhang (zhangshuang_jz@sina.cn).
