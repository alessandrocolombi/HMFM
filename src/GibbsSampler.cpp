//
// Created by pietr on 12/12/2021.
//
#include "FullConditional.h"
#include "FC_tau.h"
#include "FC_U.h"
#include "FC_S.h"
#include "FC_Lambda.h"
#include "FC_Mstar.h"
#include "FC_gamma.h"
#include "GibbsSampler.h"

GibbsSampler::GibbsSampler(unsigned int n, unsigned int b, unsigned int t, GS_data g) {
        n_iter=n;
        burn_in=b;
        thin=t;
        gs_data=g;
        FC_tau* tau;
        FC_U* U;
        FC_S* S;
        FC_Mstar* Mstar;
        FC_gamma* gamma;
        Partition* Partition;
        FC_Lambda* lambda;

    std::vector<FullConditional*> fc=FullConditionals ={Partition, tau, U, S, Mstar, gamma,lambda};
        FullConditionals=fc;
        //out={{"M*", vec}, {"K", vec}, {"U", vec}, {"S", vec},{"tau", vec},{"gamma", vec},{"adaptvarpopgamma", vec}};

}

void GibbsSampler::GS_Step() {
  for(FullConditional* full_cond: this->FullConditionals){ // mettere prima update della partition (da aggiungere anche prima)
    std::cout<<1<<std::endl;
    full_cond->update(gs_data, random_engine);
    std::cout<<full_cond->name<<std::endl;
  }
  //std::cout << gs_data.M << "\n";
  //std::cout << gs_data.K << "\n";
}


//in questo se si usa la struct perdiamo l'eleganza di questo ciclo ma al
void GibbsSampler::store_params_values() {

    out.K.push_back(gs_data.K);
    out.Mstar.push_back(gs_data.Mstar);
    out.lambda.push_back(gs_data.lambda);
    out.Ctilde.push_back(gs_data.Ctilde);
    out.S.push_back(gs_data.S);
    //out.tau.push_back(gs_data.tau);
    out.U.push_back(gs_data.U);
    out.gamma.push_back(gs_data.gamma);

}


out_data GibbsSampler::sample() {
    for(unsigned int it; it<burn_in + n_iter * thin; it++){
        this->GS_Step();
        if(it>burn_in && it%thin == 0){

            this->store_params_values();
        }
    }
    return out;
}
