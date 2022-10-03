# RPCAtree
An adaptive robust PCA algorithm for embedded tree structure recovery from noisy and heterogeneous data

## Usage

To use RPCAtree, please set `matlab` working directory to RPCAtree after cloning this repository. Then open source code `carryout_RPCAtree.m` to execute the commands line by line as following.

* Input of RPCAtree

  The input data is single-cell RNA sequencing data or single-cell DNA sequencing data (SNV data or CNV data), where the rows of input matrix represent data points, such as cells, and the columns represent features, such as genes, mutations, etc. 

```
load tree_300.mat
% Normalization
XX = (X-min(min(X)))./(max(max(X))-min(min(X)));
D = XX;
[m,n]=size(D);

%%  If there are no missing values
omega=1:(m*n);
omegaC=[];

%% If there exists missing entries
% miss=2;
% omega=find(D~=miss);
% omegaC=find(D==miss);
```
  
* Set adaptive parameter

```
lambda=1/sqrt(max(m,n))*(1+3*length(omegaC)/(m*n));
varD=var(D);
sigma= (sum(varD))/max(m,n);
gamma=m/(n*sqrt(n));
theta=sqrt(m)/(n*sqrt(n));
```

* Run RPCAtree
```
tol = 1e-10;
maxIter = 1000;
tic
[A1, E1, C, R, B] = RPCAtree(D, omega, lambda, theta, gamma, sigma, tol, maxIter);
toc
```
  

* Output of RPCAtree

  A1: the recovered low-rank matrix;

  E1: the sparse matrix;

  C: the cluster center;

  R: right stochastic matrix;

  B: the adjacency matrix of the cluster centers.

