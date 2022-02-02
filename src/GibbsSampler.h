//
// Created by pietr on 12/12/2021.
//

#ifndef GDFMM_GIBBSSAMPLER_H
#define GDFMM_GIBBSSAMPLER_H

#include <gsl/gsl_rng.h>     //For random number generators
#include <gsl/gsl_randist.h> //For random variates and probability density functions
#include <gsl/gsl_cdf.h> 	 //For cumulative density functions
#include <gsl/gsl_bspline.h> //For spline operations
#include <gsl/gsl_linalg.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "FullConditional.h"
#include "GS_data.h"
#include "out_data.h"
#include "FC_tau.h"
#include "FC_U.h"
#include "FC_S.h"
#include "FC_Lambda.h"
#include "FC_Mstar.h"
#include "FC_gamma.h"
#include <chrono>

class GibbsSampler {
public:
    unsigned int n_iter;
    unsigned int burn_in;
    unsigned int thin;
    void sample();
    // Constructor for update of all GDFMM values
    GibbsSampler(Eigen::MatrixXd const &data, unsigned int n, unsigned int b_in,
             unsigned int thn, unsigned int seed,std::string P0_prior_name, Rcpp::List option) ;
    // Constructor when number of components (M) is fixed
    GibbsSampler(Eigen::MatrixXd const &data, unsigned int n, unsigned int b_in,
            unsigned int thn, unsigned int seed, std::string P0_prior_name, unsigned int M,
            Rcpp::List option);
    // Data structure for the output
    out_data out;

private:
    std::vector<std::shared_ptr<FullConditional> > FullConditionals; // vector of shared pointer to FC class
    sample::GSL_RNG random_engine; // GSL random engine to sample from random distribution
    GS_data gs_data; // data structure to store values that are updated during Gibbs Sampler
    std::string model;
    void store_params_values();
    void GS_Step();
    bool M_fixed;
};


#endif //GDFMM_GIBBSSAMPLER_H
