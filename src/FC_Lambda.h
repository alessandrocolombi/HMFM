#ifndef __FC_LAMBDA_H__
#define __FC_LAMBDA_H__

#include "FullConditional.h"

class FC_Lambda : public FullConditional{
private:
    /* data */
    double a2;
    double b2;
public:
    FC_Lambda(/* args */);
    ~FC_Lambda();
    void update(GS_data& gs_data, sample::GSL_RNG gs_engine) override;
};

#endif