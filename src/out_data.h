#ifndef __OUTDATA_H
#define __OUTDATA_H

#include "include_headers.h"
#include "recurrent_traits.h"
struct out_data{
  /* DATA */
  // single values
  std::vector<unsigned int> K; // number of allocated component ==> number of clusters
  std::vector<unsigned int> Mstar; // number of NON-allocated component
  std::vector<double> lambda; // M|lambda ~ Poi(lambda)
  //vectors
  std::vector<std::vector< std::vector<unsigned int>>> Ctilde;
  std::vector<std::vector< std::vector<double>>> S;
  std::vector<std::vector< std::vector<double>>> tau;
  //output
  std::vector<std::vector<double>> U; // auxiliary variable
  std::vector<std::vector<double>> gamma; // vector of d gamma, one for each group
  //-----------------------------------------------------------//
  /* CONSTRUCTOR */
};

#endif
