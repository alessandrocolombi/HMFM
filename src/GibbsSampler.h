//
// Created by pietr on 12/12/2021.
//

#ifndef GDFMM_GIBBSSAMPLER_H
#define GDFMM_GIBBSSAMPLER_H

#include <gsl/gsl_rng.h>     //For random number generators
#include <gsl/gsl_randist.h> //For random variates and probability density functions
#include <gsl/gsl_cdf.h> 	 //For cumulative density functions
#include <gsl/gsl_bspline.h> //For spline operations
#include <gsl/gsl_linalg.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "FullConditional.h"
#include "GS_data.h"
#include "out_data.h"
#include "FC_tau.h"
#include "FC_U.h"
#include "FC_S.h"
#include "FC_Lambda.h"
#include "FC_Mstar.h"
#include "FC_gamma.h"
#include <chrono>



typedef std::vector<double> params;
using std::string;
// BISOGNA METTERE UN IF SEED IMPOSTATO DA UTENTE
class GibbsSampler {
public:
    unsigned int n_iter;
    unsigned int burn_in;
    unsigned int thin;
    out_data sample();
    GibbsSampler(Eigen::MatrixXd const &data, unsigned int n, unsigned int b_in,
             unsigned int thn, unsigned int seed, Rcpp::List option) ;
    out_data out;

private:
    std::vector<FullConditional*> FullConditionals; //potrebbe diventare un array? Passato
    sample::GSL_RNG random_engine;
    GS_data gs_data; // Passato
     // default -> idealmente questo era output_data che però è diventato una struct

    //BOSCA -> PETER  peter qua io ho rimesso la struct ma non sono sicuro sia il metodo migliore

    //std::map<string, std::vector<double>> output_data; // Questo diventa una struct

    //BOSCA-->PETER FERRETTI MI SA CHE BISOGNA CAMBIARE UN PO' I METODI IN MODO CHE VADANO A SALVARE I PARAMETRI NELLE STRUCT
    //se ho fatto qualcosa di poco sensato stravolgi pure tutto eh
    std::map<string, double> parameters{{"M", 0.0}, {"M*", 0.0}, {"K", 0.0}};// Forse sostituito da Gs_data
    string model;
    void store_params_values();
    void GS_Step();
    unsigned int seed;
};


#endif //GDFMM_GIBBSSAMPLER_H
