#ifndef __FC_LAMBDAMARGINAL_H__
#define __FC_LAMBDAMARGINAL_H__

#include "FC_Lambda.h"

class FC_LambdaMarginal : public FC_Lambda{
public:
    FC_LambdaMarginal(std::string _na, double _a, double _b, bool _keepfixed) : FC_Lambda(_na,_a,_b,_keepfixed){};
    FC_LambdaMarginal(std::string _na, bool _keepfixed): FC_Lambda(_na,_keepfixed){};
    FC_LambdaMarginal(bool _keepfixed):FC_Lambda(_keepfixed){};
    ~FC_LambdaMarginal() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);
};

#endif
