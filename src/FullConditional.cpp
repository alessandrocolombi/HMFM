#include "FullConditional.h"

bool FullConditional::binary_decision(double p, sample::GSL_RNG engine){
    sample::runif unif;
    double u = unif(engine);

    return u<p;
}