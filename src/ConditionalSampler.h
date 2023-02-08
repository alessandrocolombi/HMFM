#ifndef GDFMM_CONDITIONALSAMPLER_H
#define GDFMM_CONDITIONALSAMPLER_H

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>      // For progress bar
#include <progress_bar.hpp>
#include <gsl/gsl_rng.h>     //For random number generators
#include <gsl/gsl_randist.h> //For random variates and probability density functions
#include <gsl/gsl_cdf.h> 	 //For cumulative density functions
#include <gsl/gsl_bspline.h> //For spline operations
#include <gsl/gsl_linalg.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "FullConditional.h"
#include "FC_tau_mv.h"
#include "FC_U.h"
#include "FC_S.h"
#include "FC_Lambda.h"
#include "FC_Mstar.h"
#include "FC_gamma.h"
#include "FC_Partition_mv.h"
#include "GS_data.h"
#include "out_data.h"
#include "Individual.h"
#include <chrono>

class ConditionalSampler {
public:

    // sampling parameters
    unsigned int n_iter;
    unsigned int burn_in;
    unsigned int thin;

    // data parameters
    unsigned int n;
    unsigned int d;
    std::vector<unsigned int> n_j;
    std::vector<std::string> ID_i;
    std::vector<unsigned int> s_i;

    // sampling function
    void sample();

    // Constructor 
    ConditionalSampler( const Rcpp::List& _data_list, 
                        unsigned int _n_it, unsigned int _b_in, unsigned int _thn,
                        unsigned int _seed, std::string _P0_prior_name, bool _FixPart, 
                        const Rcpp::List& _option);
                
    // Data structure for the output
    out_data out;
    
private:
    std::vector<std::shared_ptr<FullConditional> > FullConditionals; // vector of shared pointer to FC class
    sample::GSL_RNG random_engine; // GSL random engine to sample from random distribution
    GS_data gs_data; // data structure to store values that are updated during Gibbs Sampler
    void store_params_values();
    void GS_Step();
    bool Partition_fixed;
};


#endif 
