#ifndef GDFMM_FC_GAMMA_H
#define GDFMM_FC_GAMMA_H
#include "FullConditional.h"
#include <gsl/gsl_randist.h>


class FC_gamma: public FullConditional {
private:
    /* MEMBERS */
    double hyp1 = 0.234;
    double hyp2 = 0.7;
    double Lambda;
    double adapt_var_pop_gamma;
    int alpha;
    int beta;
    /* METHODS */
    double log_full_gamma(double x, double Lambda, int k, double M_na, GDFMM_Traits::MatUnsCol n_jk) const;
    double  sumlgamma(double x, GDFMM_Traits::MatUnsCol n_jk) const;
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ DA RIVEDERE

public:
  std::string name="gamma";
    FC_gamma(/* args */) {};
    ~FC_gamma() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;


};


#endif //GDFMM_FC_GAMMA_H
