#ifndef GDFMM_INDIVIDUAL_H
#define GDFMM_INDIVIDUAL_H

#include "include_headers.h"
#include "recurrent_traits.h"
#include <Rcpp.h>

class Individual{
public :
    Individual(const std::string& _id, unsigned int _n, double _mean, double _var ) : ID(_id), n_ji(_n), mean_ji(_mean), var_ji(_var){};
    Individual() = default;

    std::string ID; // individual ID
    unsigned int n_ji; // number of observations taken for individual i in level j
    double mean_ji; // mean of observations taken for individual i in level j
    double var_ji; // variance of observations taken for individual i in level j

    virtual ~Individual() = default;
};
#endif
