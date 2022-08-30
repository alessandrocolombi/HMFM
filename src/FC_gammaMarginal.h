#ifndef GDFMM_FC_GAMMAMARGINAL_H
#define GDFMM_FC_GAMMAMARGINAL_H
#include "FullConditional.h"
#include "FC_gamma.h"
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_gammaMarginal: public FC_gamma {
public:

    FC_gammaMarginal(std::string _na, bool _keepfixed) : FC_gamma(_na, _keepfixed){};
    FC_gammaMarginal(std::string _na, double _h1, double _h2, double _pow, double _adapt_var0, int _a, int _b, bool _keepfixed) : FC_gamma(_na,_h1,_h2,_pow,_adapt_var0,_a,_b,_keepfixed){};
    ~FC_gammaMarginal() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);


};


#endif //GDFMM_FC_GAMMAMARGINAL_H