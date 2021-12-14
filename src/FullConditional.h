#include "include_headers.h"
#include "GDFMM_types.h"
#include "GSL_wrappers.h"

class FullConditional{
    std::map< std::string, double> hyper_param;

public :
    virtual void update(GS_parameters& param) = 0;
    bool binary_decision(double p1);
};