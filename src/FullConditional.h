#ifndef GDFMM_FULLCONDITIONAL_H
#define GDFMM_FULLCONDITIONAL_H


#include "include_headers.h"
#include "recurrent_traits.h"
#include "GS_data.h"
#include "GSL_wrappers.h"

class FullConditional{
public :
  std::string name;
    virtual void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) = 0;
    void print() const;
    bool binary_decision(double p1, const sample::GSL_RNG& engine) const;
};
#endif
