//
// Created by pietr on 12/12/2021.
//
#include "GibbsSampler.h"


GibbsSampler::GibbsSampler(Eigen::MatrixXd const &data, unsigned int n, unsigned int b, unsigned int t){
        n_iter=n;
        burn_in=b;
        thin=t;
        GS_data g(data, n,b,t);
        //pensare un modo per cambiare questo
        gs_data=g;
        /*
        FC_tau tau;
        FC_tau *ptau=&tau;
        FC_U U;
        FC_U *pU=&U;
        FC_S S;
        FC_S *pS=&S;
        FC_Mstar Mstar;
        FC_Mstar *pMstar=&Mstar;
        FC_gamma gamma;
        FC_gamma *pgamma=&gamma;
        Partition partition;
        Partition *pPartition=&partition;
        FC_Lambda lambda;
        FC_Lambda *plambda=&lambda;
        */

        FC_Mstar Mstar;
        FC_gamma gamma;
        FC_tau tau;
        FC_U U;
        FC_S S;
        Partition partition;
        FC_Lambda lambda;


        std::vector<FullConditional*> fc{&partition, &Mstar, &tau, &U, &S, &gamma,&lambda};
        FullConditionals=fc;
        //out={{"M*", vec}, {"K", vec}, {"U", vec}, {"S", vec},{"tau", vec},{"gamma", vec},{"adaptvarpopgamma", vec}};

}

void GibbsSampler::GS_Step() {
int k=0;
  for(auto full_cond : FullConditionals){
    std::cout<<"questo è il nome della partition di cui facciamo l'update "<<std::endl;// mettere prima update della partition (da aggiungere anche prima)
    //SI BLOCCA QUI, NON RIESCE A PRENDERE FULL COND
    //std::cout<<full_cond->name<<std::endl;
    std::cout<<*full_cond<<std::endl;
    //full_cond->update(gs_data, random_engine);
    k=k+1;
    std::cout<<k<<std::endl;

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
