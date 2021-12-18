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
#include "stdlib.h"
#include "stdio.h"
#include "string"
#include "vector"
#include "map"
#include "FullConditional.h"

typedef std::vector<float> params;
using std::string;

class GibbsSampler {
public:
    static unsigned int n_iter;
    unsigned int burn_in;
    unsigned int thin;

    std::map<string, std::vector<float>> sample();
    GibbsSampler(*args);
    ~GibbsSampler();

private:
    std::map<string, params(n_iter)> output_data;
    std::map<string, float> parameters;
    string model;
    std::map<string, float> initial_values;
    std::vector<FullConditional*> FullConditionals;

    void store_params_values();
    void GS_Step();
    unsigned int seed;
};


#endif //GDFMM_GIBBSSAMPLER_H
