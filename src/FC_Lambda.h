#ifndef __FC_LAMBDA_H__
#define __FC_LAMBDA_H__

#include "FullConditional.h"

class FC_Lambda : public FullConditional{
protected:
    /* hyper parameters */

    // Lambda hyperparameters: Lambda \sim gamma(a2,b2)
    double a2 = 1.0; 
    double b2 = 1.0;

    // gamma|Lambda hyperparameters: gamma_j|Lambda \sim gamma(a_gamma, b_gamma*Lambda)
    double a_gamma = 0.0;
    double b_gamma = 0.0; 
public:
    FC_Lambda(std::string na, double _a, double _b, double _agamma, double _bgamma, bool _keepfixed) : FullConditional(na,_keepfixed), a2(_a), b2(_b), a_gamma(_agamma), b_gamma(_bgamma) {};
    FC_Lambda(std::string na, double a, double b, bool _keepfixed) : FullConditional(na,_keepfixed), a2(a), b2(b){};
    FC_Lambda(std::string na, bool _keepfixed): FullConditional(na,_keepfixed){};
    FC_Lambda(bool _keepfixed):FullConditional("Lambda",_keepfixed){};
    ~FC_Lambda() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
};

#endif
