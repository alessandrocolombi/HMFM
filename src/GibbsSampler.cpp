//
// Created by pietr on 12/12/2021.
//
#include "GibbsSampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>


GibbsSampler::GibbsSampler(Eigen::MatrixXd const &data, unsigned int n, unsigned int b,
             unsigned int t, unsigned int seed):random_engine(seed){
        n_iter=n;
        burn_in=b;
        thin=t;
        GS_data g(data, n_iter, burn_in, thin, random_engine);
        //pensare un modo per cambiare questo
        //gs_data=g;
        Partition partition("Partition");
        // g.p=&partition;
        FC_Mstar Mstar("Mstar");
        FC_gamma gamma("gamma");
        FC_tau tau("tau");
        FC_U U("U");
        FC_S S("S");

        FC_Lambda lambda("lambda");
        std::vector<FullConditional*> fc{&U, &partition, &Mstar, &gamma,  &S, &tau, &lambda};
        FullConditionals=fc;
        //partition.update(gs_data, random_engine);
       // Mstar.update(gs_data, random_engine);

        //tau.update(gs_data, random_engine);
       // Rcpp::Rcout<< <<std::endl;
        //std::cout<< fc[1]->name<<std::endl;
     int k=0;   //out={{"M*", vec}, {"K", vec}, {"U", vec}, {"S", vec},{"tau", vec},{"gamma", vec},{"adaptvarpopgamma", vec}};
    for(unsigned int it=0; it<burn_in + n_iter * thin; it++){
        Rcpp::Rcout<< it<<std::endl;
        for(FullConditional* full_cond: FullConditionals){
            Rcpp::Rcout<< "Update Step : " << full_cond->name <<std::endl;

            auto t_start = std::chrono::high_resolution_clock::now();
            full_cond->update(g, random_engine);
            auto t_end = std::chrono::high_resolution_clock::now();
            double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
            Rcpp::Rcout << "It took "<< elapsed_time_ms <<" msecond(s) to update "<< full_cond->name<<std::endl;

        }
        Rcpp::Rcout<<"\nValue of M : " << g.M << " - Value of K : " << g.K <<std::endl;
        //Rcpp::Rcout<< "finoa qua"<<std::endl;
        if(it>burn_in && it%thin == 0){
            out.K.push_back(g.K);
            out.Mstar.push_back(g.Mstar);
            out.lambda.push_back(g.lambda);
            out.Ctilde.push_back(g.Ctilde);
            out.S.push_back(g.S);
            std::vector< std::vector<double>> tau;
            tau.push_back(g.mu);
            tau.push_back(g.sigma);
            out.tau.push_back(tau);
            out.U.push_back(g.U);
            out.gamma.push_back(g.gamma);
        }
        k=k+1;
       // Rcpp::Rcout<< k;
    }

}

void GibbsSampler::GS_Step() {

//std::cout<<FullConditionals[0]->name<<std::endl;
//FullConditionals[0]->update(gs_data, random_engine);
   //for(FullConditional* full_cond: this->FullConditionals){


   // std::cout<<"update solo partition"<<std::endl;// mettere prima update della partition (da aggiungere anche prima)
    //SI BLOCCA QUI, NON RIESCE A PRENDERE FULL COND
    //partition.update(
    //Rcpp::Rcout<<full_cond->name<<std::endl;



    //full_cond->update(gs_data, random_engine);
    //k=k+1;
   // std::cout<<k<<std::endl;
   //}
}

  //std::cout << gs_data.M << "\n";
  //std::cout << gs_data.K << "\n";



//in questo se si usa la struct perdiamo l'eleganza di questo ciclo ma al
void GibbsSampler::store_params_values() {

    /*out.K.push_back(g.K);
    out.Mstar.push_back(g.Mstar);
    out.lambda.push_back(gs_data.lambda);
    out.Ctilde.push_back(gs_data.Ctilde);
    out.S.push_back(gs_data.S);
    //out.tau.push_back(gs_data.tau);
    out.U.push_back(gs_data.U);
    out.gamma.push_back(gs_data.gamma);*/

}


out_data GibbsSampler::sample() {
    for(unsigned int it=0; it<burn_in + n_iter * thin; it++){

        this->GS_Step();
        if(it>burn_in && it%thin == 0){

            this->store_params_values();
        }
    }
    return out;
}
