
#ifndef GDFMM_FC_GAMMA_H
#define GDFMM_FC_GAMMA_H
#include "FullConditional.h"
#include <gsl/gsl_randist.h>


class FC_gamma: public FullConditional {
private:
    float hyp1=0.234;
    float hyp2=0.7;
    double Mna = Gs_data::Mstar;
    double iter=GS_data::iterations;
    double Lambda;
    unsigned int K=GS_data::K;
    float adapt_pop_gamma;
    int alpha;
    int beta;
    double  log_full_gamma(double x, double Lambda, int k, double M_na, double n_jk);
    int a1;
    int b1;

public:
    FC_gamma(/* args */);
    ~FC_gamma();
    void update(GS_data& gs_data, sample::GSL_RNG gs_engine) override;


};


#endif //GDFMM_FC_GAMMA_H
