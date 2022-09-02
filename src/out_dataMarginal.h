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
};

#endif
