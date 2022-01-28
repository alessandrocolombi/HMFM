#ifndef GDFMM_FC_GAMMA_H
#define GDFMM_FC_GAMMA_H
#include "FullConditional.h"
#include <gsl/gsl_randist.h>
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_gamma: public FullConditional {
private:
    /* MEMBERS */
    double hyp1 = 0.234;
    double hyp2 = 0.7;
    unsigned int power = 2;
    double adapt_var_pop_gamma=1;//valore iniziale che viene aggiornato
    int alpha = 1;
    int beta = 1;
    /* METHODS */
    double log_full_gamma(double x, double Lambda, unsigned int k, unsigned int M_star, const GDFMM_Traits::MatUnsCol & n_jk);
    double sumlgamma(double x, const GDFMM_Traits::MatUnsCol& n_jk);
    double l_dgamma(double gamma, double a, double b);

public:

    FC_gamma(std::string na){name=na;};
    FC_gamma(std::string na, double h1, double h2, double pow, double adapt_var0, int a, int b) : hyp1(h1),
            hyp2(h2), power(pow), adapt_var_pop_gamma(adapt_var0), alpha(a), beta(b) {name = na;};
    ~FC_gamma() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;


};


#endif //GDFMM_FC_GAMMA_H
