#include "FullConditional.h"
bool binary_decision(double p1, const sample::GSL_RNG& engine) const;
bool FullConditional::binary_decision(double p1, const sample::GSL_RNG& engine) const{
    // Sample u from a Uniform(0,1) and verify if is less than p
    sample::runif unif;
    double u = unif(engine);

    return u<p;
}