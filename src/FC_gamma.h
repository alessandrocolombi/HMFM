
#ifndef GDFMM_FC_GAMMA_H
#define GDFMM_FC_GAMMA_H
#include "FullConditional.h"
#include <gsl/gsl_randist.h>


class FC_gamma: public FullConditional {
private:
    /* MEMBERS */
    double hyp1 = 0.234;
    double hyp2 = 0.7;
    double Mna = Gs_data::Mstar;
    double iter = GS_data::iterations;
    double Lambda;
    unsigned int K = GS_data::K;
    double adapt_var_pop_gamma;
    int alpha;
    int beta;
    /* METHODS */
    double log_full_gamma(double x, double Lambda, int k, double M_na, double n_jk) const; 
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ DA RIVEDERE
        
public:
    FC_gamma(/* args */);
    ~FC_gamma();
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;


};


#endif //GDFMM_FC_GAMMA_H
