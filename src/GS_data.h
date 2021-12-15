#ifndef __GS_DATA_H
#define __GS_DATA_H

#include "include_headers.h"

struct GS_data{
    // single values
    unsigned int d; // number of groups
    unsigned int k; // number of allocated component ==> number of clusters
    unsigned int Mstar; // number of NON-allocated component
    unsigned int M; // total number of component
    double lambda; // M|lambda ~ Poi(lambda)
    // vectors
    std::vector<double> U; // auxiliary variable
    std::vector<double> gamma; // vector of d gamma, one for eaach group
    // matrix or vector of vectors
    
};

#endif