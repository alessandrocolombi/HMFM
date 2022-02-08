#ifndef __OUTDATA_H
#define __OUTDATA_H

#include "include_headers.h"
#include "recurrent_traits.h"

struct out_data{
  /* DATA */
  // PARTITION
  std::vector<unsigned int> K; // number of allocated component ==> number of clusters
  std::vector<unsigned int> Mstar; // number of NON-allocated component
  std::vector<std::vector< std::vector<unsigned int>>> Ctilde;

  // PARAMETERS OF MIXTURES
  std::vector< std::vector<double>> mu;
  std::vector< std::vector<double>> sigma;
  std::vector<GDFMM_Traits::MatRow> w_ji;

  // ALGORITHIMC INTERESTING VALUES
  std::vector<double> lambda; // M|lambda ~ Poi(lambda)
  std::vector<double> U; // auxiliary variable
  std::vector<double> gamma; // vector of d gamma, one for each group
  std::vector<double> adapt_var;
  //-----------------------------------------------------------//
  /* CONSTRUCTOR */
};

#endif
