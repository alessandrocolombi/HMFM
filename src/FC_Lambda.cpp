#include "FC_Lambda.h"

void FC_Lambda::update(GS_parameters& param, sample::GSL_RNG engine){
    // From param all needed variable are retrived
    int k = param["k"][0][0];
    int d = param["c"].size();
    double gamma = param["gamma"][0]; // CHIEDERE A COLMBI SE ABBIAMO d gamma UGUALI O d gamma[j]
    std::vector<double> U = param["U"][0];
    // Random sampler are created
    sample::runif Unif;
    sample::rgamma Gamma;
    // Update routine
    int a2_star = d*(k-1) + a2; // a2 può essere un double?-->nel caso dovrei trattare questa operazione

    double sum = 0.0;

    for(const double & u : U){
        sum += // RIPRENDI DA QUI   <==================================================================
    }

}