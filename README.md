# RPCAtree
An adaptive robust PCA algorithm for embedded tree structure recovery from noisy and heterogeneous data

## Usage
* Input of RPCAtree

  The input data is single-cell RNA sequencing data or single-cell DNA sequencing data (SNV data or CNV data), where the rows of input matrix represent data points, such as cells, and the columns represent features, such as genes, mutations, etc. 

* Output of RPCAtree

  A1: the recovered low-rank matrix;

  E1: the sparse matrix;

  C: the cluster center;

  R: right stochastic matrix;

  B: the adjacency matrix of the cluster centers.

