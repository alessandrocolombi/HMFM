#include "include_headers.h"
#include "GS_data.h"
#include "GSL_wrappers.h"

class FullConditional{
    
public :
    virtual void update(GS_data& gs_data, sample::GSL_RNG gs_engine) = 0;
    bool binary_decision(double p1, sample::GSL_RNG engine);
};