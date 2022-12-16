#ifndef GDFMM_FC_GAMMAMARGINAL_H
#define GDFMM_FC_GAMMAMARGINAL_H
#include "FullConditional.h"
#include "FC_gamma.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <gsl/gsl_sf_psi.h> // for digamma function

class FC_gammaMarginal: public FC_gamma {
private:
    // evaluate log_pi(gamma_1,...,gamma_d = x_1,...,x_d | rest )
    double log_FCgamma_marginal(const std::vector<double>& x, const double& Lambda, const unsigned int& K,
                                const std::vector<double>& U, const GDFMM_Traits::MatUnsCol& N) const; 
    // evaluate grad_log_pi(gamma_1,...,gamma_d = x_1,...,x_d | rest )
    std::vector<double> grad_log_FCgamma_marginal(const std::vector<double>& x, const double& Lambda, const unsigned int& K, 
                                                  const std::vector<double>& U, const GDFMM_Traits::MatUnsCol& N) const; 
    // compute log( (a)_n )   
    double log_raising_factorial(const unsigned int& n, const double& a)const;
public:
    FC_gammaMarginal(std::string _na, double _h1, double _h2, double _pow, unsigned int _d, double _adapt_var0, int _a, int _b, bool _keepfixed) : 
                        FC_gamma(_na,_h1,_h2,_pow,_d,_adapt_var0,_a,_b,_keepfixed){};
    ~FC_gammaMarginal() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);


};


#endif //GDFMM_FC_GAMMAMARGINAL_H