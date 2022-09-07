#ifndef __FC_LAMBDAMARGINAL_H__
#define __FC_LAMBDAMARGINAL_H__

#include "FC_Lambda.h"

class FC_LambdaMarginal : public FC_Lambda{
private:
    double hyp1 = 0.234;
    double hyp2 = 0.7;
    unsigned int power = 10;
    double adapt_var_proposal_Lambda; //variance to be adapted in MH steps
    double log_FCLambda_marginal(const double& x, const std::vector<double>& U, const std::vector<double>& gamma, const unsigned int& K) const;
public:
    FC_LambdaMarginal(std::string _na, double _a, double _b, bool _keepfixed, double _h1, double _h2, double _pow, double _adapt_var0);
    //FC_LambdaMarginal(std::string _na, bool _keepfixed): FC_Lambda(_na,_keepfixed){};
    //FC_LambdaMarginal(bool _keepfixed):FC_Lambda(_keepfixed){};
    ~FC_LambdaMarginal() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);
};

#endif
