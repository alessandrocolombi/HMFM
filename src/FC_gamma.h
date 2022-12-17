#ifndef GDFMM_FC_GAMMA_H
#define GDFMM_FC_GAMMA_H
#include "FullConditional.h"
#include <gsl/gsl_randist.h>
#include <Rcpp.h>
#include <RcppEigen.h>



class FC_gamma: public FullConditional {
protected:
    /* MEMBERS */
    double hyp1 = 0.234;
    double hyp2 = 0.7;
    double s_p  = 0.01; // mala parameter
    unsigned int power = 10;
    std::vector<double> adapt_var_pop_gamma; //vector for variances to be adapted in MH steps
    int alpha = 1;
    int beta = 1;
    /* METHODS */
    double log_full_gamma(double x, double Lambda, unsigned int k, unsigned int M_star, const GDFMM_Traits::MatUnsCol & n_jk);
    double sumlgamma(double x, const GDFMM_Traits::MatUnsCol& n_jk);
    double l_dgamma(double gamma, double a, double b);

public:

    //FC_gamma(std::string na, bool _keepfixed) : FullConditional(na,_keepfixed){};
    FC_gamma(std::string na, double h1, double h2, double pow, unsigned int d, double adapt_var0, int a, int b, double _s_p, bool _keepfixed);
    FC_gamma(std::string na, double h1, double h2, double pow, unsigned int d, double adapt_var0, int a, int b, bool _keepfixed):
                    FC_gamma(na,h1,h2,pow,d,adapt_var0,a,b,0.01,_keepfixed){};
    ~FC_gamma() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;


};


#endif //GDFMM_FC_GAMMA_H
