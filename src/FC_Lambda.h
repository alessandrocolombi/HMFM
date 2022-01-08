#ifndef __FC_LAMBDA_H__
#define __FC_LAMBDA_H__

#include "FullConditional.h"

class FC_Lambda : public FullConditional{
private:
    /* hyper parameters */
    double a2;
    double b2;
public:

    FC_Lambda(/* args */) {};
    FC_Lambda(std::string na){name=na;};
    FC_Lambda(double a, double b): a2(a), b2(b){};
    ~FC_Lambda() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
};

#endif
