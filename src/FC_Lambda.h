#ifndef __FC_LAMBDA_H__
#define __FC_LAMBDA_H__

#include "FullConditional.h"

class FC_Lambda : public FullConditional{
private:
    /* hyper parameters */
    double a2 = 1;
    double b2 = 1;
public:
    FC_Lambda(std::string na, double a, double b, bool _keepfixed) : FullConditional(na,_keepfixed), a2(a), b2(b){}//, keep_fixed(_keepfixed){name = na;};
    FC_Lambda(std::string na, bool _keepfixed): FullConditional(na,_keepfixed){};
    FC_Lambda(bool _keepfixed):FullConditional("Lambda",_keepfixed){};//keep_fixed(_keepfixed){name = "Lambda";};
    ~FC_Lambda() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
};

#endif
