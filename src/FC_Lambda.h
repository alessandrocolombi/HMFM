#ifndef __FC_LAMBDA_H__
#define __FC_LAMBDA_H__

#include "FullConditional.h"

class FC_Lambda : public FullConditional{
private:
    /* data */
    a2;
    b2;
public:
    FC_Lambda(/* args */);
    ~FC_Lambda();
    void update(GS_parameters& param, sample::GSL_RNG engine) override;
};

#endif