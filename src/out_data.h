#ifndef __OUTDATA_H
#define __OUTDATA_H

#include "include_headers.h"
#include "recurrent_traits.h"
#include <Rcpp.h>
#include <RcppEigen.h>

struct out_data{
  /* DATA */
  // PARTITION
  std::vector<unsigned int> K; // number of allocated component ==> number of clusters, length is equal to n_iter
  std::vector<unsigned int> Mstar; // number of NON-allocated component
  std::vector<std::vector< std::vector<unsigned int>>> Ctilde;

  // PARAMETERS OF MIXTURES
  /* Old version - ragazzi
  std::vector< std::vector<double>> mu; 
  std::vector< std::vector<double>> sigma;
  */
  std::vector< Rcpp::NumericVector > mu; 
  std::vector< Rcpp::NumericVector > sigma;
  std::vector< GDFMM_Traits::MatRow > S;
  std::vector< GDFMM_Traits::MatRow > beta;
  // if partition is fixed estimated weights for the components are stored in w_jk
  std::vector<GDFMM_Traits::MatRow> w_jk; //POSSO TOGLIERE COMPLETAMENTE IL CALCOLO DEI w_jk??

  // ALGORITHIMC INTERESTING VALUES
  std::vector<double> lambda; // M|lambda ~ Poi(lambda)
  std::vector<double> U; // auxiliary variable
  std::vector<double> gamma; // vector of d gamma, one for each group
  std::vector<double> adapt_var;
  std::vector<double> log_prod_psiU;
  //-----------------------------------------------------------//
  /* CONSTRUCTOR */
};

#endif
