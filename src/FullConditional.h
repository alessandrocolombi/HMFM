#include "include_headers.h"
#include "GDFMM_types.h"
#include "GSL_wrappers.h"

class FullConditional{
    
public :
    virtual void update(GS_parameters& param, sample::GSL_RNG engine) = 0;
    bool binary_decision(double p1, sample::GSL_RNG engine);
};