//
// Created by pietr on 12/12/2021.
//
#include "GibbsSampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>

// Consrtructor
GibbsSampler::GibbsSampler(Eigen::MatrixXd const &data, unsigned int n_it, unsigned int b_in,
            unsigned int thn, unsigned int seed, std::string P0_prior_name, bool FixPart,
            Rcpp::List option) : random_engine(seed), Partition_fixed(FixPart){
    // Extract hyper_parameter and initialization values from option
    n_iter = n_it;
    burn_in = b_in;
    thin = thn;

    if(P0_prior_name == "Normal-InvGamma"){
        
        // Read all hyper-parameters passed with option
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
        
        std::vector<unsigned int> partition_vec;
        if(Partition_fixed){
            partition_vec = Rcpp::as<std::vector<unsigned int>>(option["partition"]);
        }
        
        // Initialize gs_data with correct random seed, given Mstar and all data assigned to same cluster
        gs_data = GS_data(data, n_iter, burn_in, thin, random_engine,
                        Mstar0, Lambda0, mu0, nu0, sigma0, P0_prior_name, partition_vec);
            
        
        //Initialize Full Conditional Objects
        auto Partition_ptr = std::make_shared<Partition>("Partition", gs_data.d, gs_data.n_j, FixPart);
        auto Mstar_ptr = std::make_shared<FC_Mstar>("Mstar", FixPart);
        auto gamma_ptr = std::make_shared<FC_gamma>("gamma", h1, h2, pow, adapt_var0, a1, b1);
        auto tau_ptr = std::make_shared<FC_tau>("tau", nu0, sigma0, mu0, k0);
        auto U_ptr = std::make_shared<FC_U>("U");
        auto S_ptr = std::make_shared<FC_S>("S");
        auto lambda_ptr = std::make_shared<FC_Lambda>("lambda", a2, b2);

        //Full Conditional vector that we will loop
        std::vector< std::shared_ptr<FullConditional> > fc{U_ptr,
                                                            Partition_ptr,
                                                            Mstar_ptr,
                                                            gamma_ptr,
                                                            S_ptr,
                                                            tau_ptr,
                                                            lambda_ptr
                                                            };
    
        std::swap(FullConditionals, fc);
    }
}

void GibbsSampler::sample() {
    for(unsigned int it = 0; it <= burn_in + n_iter*thin; it++){
        this->GS_Step();
        if(it>burn_in && it%thin == 0){
            this->store_params_values();
        }
        //updating number of iterations necessary for MH algorithm
        gs_data.iterations = it;
        if(!Partition_fixed){
            Rcpp::Rcout<< "\nIn this iteration we obtain K: "<< gs_data.K << " M: " << gs_data.M
                        <<"\n"<<std::endl;
        }
    }
}

void GibbsSampler::GS_Step() {
    //Loop for updating every fullconditional
    for(auto full_cond: FullConditionals){
        Rcpp::Rcout<< "Update Step : " << full_cond->name <<std::endl;
        //starting timer to measure updating time
        auto t_start = std::chrono::high_resolution_clock::now();
        full_cond->update(gs_data, random_engine);
        //ending timer to measure updating time
        auto t_end = std::chrono::high_resolution_clock::now();
        //elapsed time in ms
        double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        Rcpp::Rcout << "It took "<< elapsed_time_ms <<" msecond(s) to update "<< full_cond->name<<std::endl;
    }
}

void GibbsSampler::store_params_values() {
    
    if(!Partition_fixed){
        out.K.push_back(gs_data.K);
        out.Mstar.push_back(gs_data.Mstar);
        out.Ctilde.push_back(gs_data.Ctilde);
    }

    out.lambda.push_back(gs_data.lambda);
    out.S.push_back(gs_data.S);
    std::vector< std::vector<double>> tau;
    tau.push_back(gs_data.mu);
    tau.push_back(gs_data.sigma);
    out.tau.push_back(tau);
    out.U.push_back(gs_data.U);
    out.gamma.push_back(gs_data.gamma);
}
