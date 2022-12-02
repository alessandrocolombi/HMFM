#ifndef __OUTDATAMARGINAL_H
#define __OUTDATAMARGINAL_H

#include "include_headers.h"
#include "recurrent_traits.h"
#include <Rcpp.h>
#include <RcppEigen.h>

struct out_dataMarginal{

  // COUNTER
  unsigned int it_saved{0};
  // PARTITION
  std::vector<unsigned int> K; // number of allocated component ==> number of clusters
  
  //Partition[it,] contains the cluster assignments for each data in each level. Values are inserted level by level, 
  //hence Partition[it,1:n_1] are the assignments for level 1 (n_1 is the number of observations in level 1), 
  // Partition[it,(n_1+1):(n_1+n_2)] the assignments for level 2 and so on
  Rcpp::NumericMatrix Partition; // matrix of size n_iter x n_data, where n_data is the total number of data. 

  // PARAMETERS OF MIXTURES
  std::vector< Rcpp::NumericVector > mu; 
  std::vector< Rcpp::NumericVector > sigma;

  // ALGORITHIMC INTERESTING VALUES
  std::vector<double> lambda; // M|lambda ~ Poi(lambda)
  Rcpp::NumericMatrix U;      // Rcpp matrix of size d x n_iter
  Rcpp::NumericMatrix gamma;  // Rcpp matrix of size d x n_iter

  // WEIGHTS FOR DENSITY ESTIMATION
  // q is a vector of matrices. The external vector has length n_iter. Then, q[it] is a matrix of size d x K+1, where 
  // K = K[it] is the number of clusters during iteration it. Within each row, elements are the unnormalized cluster 
  // probabilities in log scale for each level j, hence 
  // log( q_1(n_11,...,n_1K,U_1) ), ... , log( q_K(n_11,...,n_1K,U_1) ), log( q_{K+1}(n_11,...,n_1K,U_1) )
  // log( q_1(n_j1,...,n_jK,U_j) ), ... , log( q_K(n_j1,...,n_jK,U_j) ), log( q_{K+1}(n_j1,...,n_jK,U_j) )
  // log( q_1(n_d1,...,n_dK,U_d) ), ... , log( q_K(n_d1,...,n_dK,U_d) ), log( q_{K+1}(n_d1,...,n_dK,U_d) )
  std::vector< GDFMM_Traits::MatRow > log_q;
};

#endif
