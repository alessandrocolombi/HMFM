//
// Created by pietr on 12/12/2021.
//
#include "GibbsSampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>


GibbsSampler::GibbsSampler(Eigen::MatrixXd const &data, unsigned int n, unsigned int b_in,
             unsigned int thn, unsigned int seed, Rcpp::List option) : random_engine(seed) {
        // Extract hyper_parameter and initialization values from option
        n_iter = n;
        burn_in = b_in;
        thin = thn;

        unsigned int Mstar0 = Rcpp::as<unsigned int>(option["Mstar0"]);
        double Lambda0 = Rcpp::as<double>(option["Lambda0"]);
        double mu0 = Rcpp::as<double>(option["mu0"]);
        double nu0 = Rcpp::as<double>(option["nu0"]);
        double sigma0 = Rcpp::as<double>(option["sigma0"]);
        double h1 = Rcpp::as<double>(option["Adapt_MH_hyp1"]);
        double h2 = Rcpp::as<double>(option["Adapt_MH_hyp2"]);
        unsigned int pow = Rcpp::as<unsigned int>(option["Adapt_MH_power_lim"]);
        double adapt_var0 = Rcpp::as<double>(option["Adapt_MH_var0"]);
        double k0 = Rcpp::as<double>(option["k0"]);
        double a1 = Rcpp::as<double>(option["alpha_gamma"]);
        double b1 = Rcpp::as<double>(option["beta_gamma"]);
        double a2 = Rcpp::as<double>(option["alpha_lambda"]);
        double b2 = Rcpp::as<double>(option["beta_lambda"]);
        
        // Initialize gs_data with the correct random seed
        gs_data = GS_data(data, n_iter, burn_in, thin, random_engine, Mstar0, Lambda0, mu0, nu0, sigma0);
        Partition partition("Partition");
        FC_Mstar Mstar("Mstar");
        FC_gamma gamma("gamma", h1, h2, pow, adapt_var0, a1, b1);
        FC_tau tau("tau", nu0, sigma0, mu0, k0);
        FC_U U("U");
        FC_S S("S");
        FC_Lambda lambda("lambda", a2, b2);

        std::vector<FullConditional*> fc{&U,
                                         &partition,
                                         &Mstar,
                                         &gamma,
                                         &S,
                                         &tau,
                                         &lambda
                                         };
        FullConditionals = fc;
        //partition.update(gs_data, random_engine);
       // Mstar.update(gs_data, random_engine);

        //tau.update(gs_data, random_engine);
       // Rcpp::Rcout<< <<std::endl;
        //std::cout<< fc[1]->name<<std::endl;
       //out={{"M*", vec}, {"K", vec}, {"U", vec}, {"S", vec},{"tau", vec},{"gamma", vec},{"adaptvarpopgamma", vec}};
    for(unsigned int it=0; it<burn_in + n_iter * thin; it++){
       // Rcpp::Rcout<< it<<std::endl;
        for(FullConditional* full_cond: FullConditionals){
            Rcpp::Rcout<< "Update Step : " << full_cond->name <<std::endl;

            auto t_start = std::chrono::high_resolution_clock::now();
            full_cond->update(gs_data, random_engine);
            auto t_end = std::chrono::high_resolution_clock::now();
            double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
            Rcpp::Rcout << "It took "<< elapsed_time_ms <<" msecond(s) to update "<< full_cond->name<<std::endl;
            gs_data.iterations = it;



        }
        Rcpp::Rcout<<"\nValue of M : " << gs_data.M << " - Value of K : " << gs_data.K << "\n"<<std::endl;
        //Rcpp::Rcout<< "finoa qua"<<std::endl;
        if(it>burn_in && it%thin == 0){
            out.K.push_back(gs_data.K);
            out.Mstar.push_back(gs_data.Mstar);
            out.lambda.push_back(gs_data.lambda);
            out.Ctilde.push_back(gs_data.Ctilde);
            out.S.push_back(gs_data.S);
            std::vector< std::vector<double>> tau;
            tau.push_back(gs_data.mu);
            tau.push_back(gs_data.sigma);
            out.tau.push_back(tau);
            out.U.push_back(gs_data.U);
            out.gamma.push_back(gs_data.gamma);
        }

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
