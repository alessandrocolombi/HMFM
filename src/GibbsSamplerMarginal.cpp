//
// Created by pietr on 12/12/2021.
//
#include "GibbsSamplerMarginal.h"
#include <Rcpp.h>
#include <RcppEigen.h>

// Consrtructor
GibbsSamplerMarginal::GibbsSamplerMarginal( Eigen::MatrixXd const &data, unsigned int n_it, unsigned int b_in,
                                            unsigned int thn, unsigned int seed, std::string P0_prior_name, bool FixPart,
                                            Rcpp::List option) : random_engine(seed), Partition_fixed(FixPart){
    // Extract hyper_parameter and initialization values from option
    n_iter = n_it;
    burn_in = b_in;
    thin = thn;

    if(P0_prior_name == "Normal-InvGamma"){

        Rcpp::Rcout<<"GibbsSamplerMarginal.cpp - Initialization"<<std::endl;

        // Manage cases if Partition is fixed
        std::vector<unsigned int> partition_vec{ Rcpp::as<std::vector<unsigned int>>(option["partition"]) }; 
        
        // Stampe
        /*
        Rcpp::Rcout<<"Partition: "<<std::endl;
        for(unsigned int i = 0; i<partition_vec.size(); i++)
            Rcpp::Rcout<<partition_vec[i]<<"; ";
        
        Rcpp::Rcout<<std::endl;
        */

        // To avoid rewriting a lot of code, I can set Mstar = 0. If so, M = K
        double Mstar0{0.0};
        //double nu{1.0}; // nu must be eliminated

        // Read all hyper-parameters passed with option
        double Lambda0 = Rcpp::as<double>(option["Lambda0"]);
        double mu0 = Rcpp::as<double>(option["mu0"]);
        double nu0 = Rcpp::as<double>(option["nu0"]);
        double sigma0 = Rcpp::as<double>(option["sigma0"]);
        double gamma0 = Rcpp::as<double>(option["gamma0"]);
        double h1 = Rcpp::as<double>(option["Adapt_MH_hyp1"]);
        double h2 = Rcpp::as<double>(option["Adapt_MH_hyp2"]);
        double s_p_U = Rcpp::as<double>(option["sp_mala_U"]);
        double s_p_gamma = Rcpp::as<double>(option["sp_mala_gamma"]);
        double k0 = Rcpp::as<double>(option["k0"]);
        double a1 = Rcpp::as<double>(option["alpha_gamma"]);
        double b1 = Rcpp::as<double>(option["beta_gamma"]);
        double a2 = Rcpp::as<double>(option["alpha_lambda"]);
        double b2 = Rcpp::as<double>(option["beta_lambda"]);
        std::vector<double> init_mean_clus{ Rcpp::as<std::vector<double>>(option["init_mean_cluster"]) };
        std::vector<double> init_var_clus{ Rcpp::as<std::vector<double>>(option["init_var_cluster"]) };
        bool FixedU = !Rcpp::as<bool>(option["UpdateU"]);
        bool FixedGamma = !Rcpp::as<bool>(option["UpdateGamma"]);
        bool FixedTau = !Rcpp::as<bool>(option["UpdateTau"]);
        bool FixedLambda = !Rcpp::as<bool>(option["UpdateLambda"]);

        // Initialize gs_data with correct random seed, given Mstar and all data assigned to same cluster
        gs_data = GS_data( data, n_iter, burn_in, thin, random_engine,
                           Mstar0, Lambda0, mu0, nu0, sigma0, gamma0, 
                           init_mean_clus, init_var_clus,
                           P0_prior_name, partition_vec);

        // Iinitialize sums in clusters
        gs_data.initialize_sums_in_clusters();
        gs_data.compute_log_prob_marginal_data(nu0, sigma0, mu0,  k0);

        //Initialize Full Conditional Objects
        auto PartitionMarginal_ptr = std::make_shared<FC_PartitionMarginal>("Partition", gs_data.d, gs_data.n_j, FixPart, nu0, sigma0, mu0, k0);
        auto gammaMarginal_ptr = std::make_shared<FC_gammaMarginal>("gamma", h1, h2, 10, gs_data.d, 1.0, a1, b1, s_p_gamma, FixedGamma);
        auto tau_ptr = std::make_shared<FC_tau>("tau", nu0, sigma0, mu0, k0, FixedTau);
        auto UMarginal_ptr = std::make_shared<FC_UMarginal>("U", FixedU, h1, h2, 10, gs_data.d, 1.0,s_p_U);
                //auto lambdaMarginal_ptr = std::make_shared<FC_LambdaMarginal>("lambda", a2, b2, FixedLambda, h1, h2, 10, 1.0); //update of lambda in conditional and marginal sampler must be the same
        auto lambda_ptr = std::make_shared<FC_Lambda>("lambda", a2, b2, FixedLambda);

        //Full Conditional vector that we will loop
        std::vector< std::shared_ptr<FullConditional> > fc{ lambda_ptr,
                                                            gammaMarginal_ptr,
                                                            UMarginal_ptr,
                                                            PartitionMarginal_ptr,
                                                            tau_ptr // tau must be update after the partition, otherwise saved values may be misleading
                                                            };

        std::swap(FullConditionals, fc);

        // Initialize return structures 
        out.K.resize(n_iter);
        out.mu.resize(n_iter);
        out.sigma.resize(n_iter);
        out.lambda.resize(n_iter);
        out.U = Rcpp::NumericMatrix(gs_data.d,n_iter);
        out.gamma = Rcpp::NumericMatrix(gs_data.d,n_iter);
        out.Partition = Rcpp::NumericMatrix(n_iter, std::accumulate(gs_data.n_j.cbegin(), gs_data.n_j.cend(), 0) ); // here; i am also computing the total number of data by summing the elements of n_j
        out.log_q.reserve(n_iter);
        out.it_saved = 0;
    }
    else{
        throw std::runtime_error("Error, P0_prior_name must be equal to Normal-InvGamma. No other cases have been implemented.");
    }

    //Rcpp::Rcout<<"DEVO CAMBIARE LA STRUTTURA DEI INPUT DEI DATI COME HA FATTO VEDERE BARTOLUCCI"<<std::endl;
}

void GibbsSamplerMarginal::sample() {

    Progress progress_bar(burn_in + n_iter*thin, TRUE); // Initialize progress bar
    for(unsigned int it = 0; it < burn_in + n_iter*thin; it++){

        // If we are in the right iteration store needed values
        // If burn_in is set to zero, the first values to be saved are the initial values.
        if(it >= burn_in && (it-burn_in)%thin == 0){
            //Rcpp::Rcout<<"it = "<<it<<std::endl;
            
            this->store_params_values();
            out.it_saved++;
        }

        // Sample from all full conditional
        this->GS_Step();

        //updating number of iterations necessary for MH algorithm
        gs_data.iterations = it;

        progress_bar.increment(); //update progress bar
    }
}

void GibbsSamplerMarginal::GS_Step() {
    //Loop for updating every fullconditional
    //Rcpp::Rcout<<"Mstar = "<<gs_data.Mstar<<"; K = "<<gs_data.K<<std::endl;
    for(auto full_cond: FullConditionals){
         
         //Rcpp::Rcout<<"-------------------------------------------"<<std::endl;
         //Rcpp::Rcout<<"gs_data.Mstar:"<<std::endl<<gs_data.Mstar<<std::endl;
         //Rcpp::Rcout<<"gs_data.K:"<<std::endl<<gs_data.K<<std::endl;
         //Rcpp::Rcout<<"gs_data.M:"<<std::endl<<gs_data.M<<std::endl;
         //Rcpp::Rcout<< "Update Step : " << full_cond->name <<std::endl;

        // starting timer to measure updating time
        // auto t_start = std::chrono::high_resolution_clock::now();
        if(!full_cond->keep_fixed){
            full_cond->update(gs_data, random_engine);
            //Rcpp::Rcout<<" --> done! "<<std::endl;
        }
        else if(full_cond->name.compare("Mstar") == 0){ //if they are equal
            throw std::runtime_error("Error, Marginal sampler must not update Mstar. ");
        }

        //Rcpp::Rcout<<"gs_data.Mstar:"<<std::endl<<gs_data.Mstar<<std::endl;
        //Rcpp::Rcout<<"gs_data.K:"<<std::endl<<gs_data.K<<std::endl;
        //Rcpp::Rcout<<"gs_data.M:"<<std::endl<<gs_data.M<<std::endl;
        //Rcpp::Rcout<<"********************************************"<<std::endl;

        // ending timer to measure updating time
        // auto t_end = std::chrono::high_resolution_clock::now();
        // elapsed time in ms
        // double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
        // Rcpp::Rcout << "It took "<< elapsed_time_ms <<" msecond(s) to update "<< full_cond->name<<std::endl;

        //Check for User Interruption
        try{
            Rcpp::checkUserInterrupt();
        }
        catch(Rcpp::internal::InterruptedException e){ 
            //Print error and return
            throw std::runtime_error("Execution stopped by the user");
        }
    }
}


void GibbsSamplerMarginal::store_params_values() {

    const unsigned int& it_saved = out.it_saved;

    if(gs_data.K==0)
        throw std::runtime_error("K is 0, this should be impossible");
    out.K[it_saved] = gs_data.K; //number of clusters

    // Save tau
    out.mu[it_saved] = Rcpp::NumericVector ( gs_data.mu.begin(),gs_data.mu.end() ); 
    out.sigma[it_saved] =   Rcpp::NumericVector (gs_data.sigma.begin(),gs_data.sigma.end()); 
    
    // Save lambda - U - gamma
    out.lambda[it_saved] = gs_data.lambda;
    out.U.column(it_saved) = Rcpp::NumericVector ( gs_data.U.begin(),gs_data.U.end() );
    out.gamma.column(it_saved) = Rcpp::NumericVector ( gs_data.gamma.begin(),gs_data.gamma.end() );
    
    // Save Partition
    std::vector<unsigned int> get_all_labels; //vector of length equal to the total number of data with all the cluster assigments at the current iteration
    for(std::size_t j = 0; j < gs_data.d; j++)
        get_all_labels.insert(get_all_labels.end(), gs_data.Ctilde[j].begin(), gs_data.Ctilde[j].end());

    out.Partition(it_saved, Rcpp::_ ) = Rcpp::NumericMatrix(    1, std::accumulate(gs_data.n_j.cbegin(), gs_data.n_j.cend(), 0), 
                                                                get_all_labels.begin() 
                                                            ); //in R notation, this is equal to Partition[it,:] = get_all_labels 

    // Save log_q
    const unsigned int& K = gs_data.K; // number of cluster
    const unsigned int& Lambda = gs_data.lambda; 
    GDFMM_Traits::MatRow Weights_it(GDFMM_Traits::MatRow::Constant(gs_data.d,gs_data.K+1,0.0)); //define d x K+1 empty matrix

    // Compute product term that is common to all levels j

    /*    
    // OLD IMPLEMENTATION THAT IS WRONG
    double log_const_new_cluster = std::inner_product(  gs_data.U.cbegin(), gs_data.U.cend(), gs_data.gamma.cbegin(), 
                                                        0.0, std::plus<>(),
                                                        [&Lambda, &K](const double& U_h, const double& gamma_h)
                                                        {  
                                                            const double Psi_h{ 1/std::pow(1.0+U_h,gamma_h) }; // compute Psi_h
                                                            return (  -gamma_h * std::log(U_h+1.0) + 
                                                                        std::log( (double)K + 1.0 + Lambda * Psi_h ) - 
                                                                        std::log( (double)K + Lambda * Psi_h)
                                                                    );
                                                        }
                                                    );
    */
    double log_const_new_cluster =  -gs_data.log_sum + 
                                    std::log( (double)K + 1.0 + Lambda*std::exp(-gs_data.log_sum) ) -
                                    std::log( (double)K + Lambda*std::exp(-gs_data.log_sum) )   ;

    // Start loop over all levels and all clusters + 1
    for(unsigned int j=0; j<gs_data.d; j++){ //for each level
        for(unsigned int l=0; l<(gs_data.K+1); l++){ // for each cluster + 1
            
            if(l<gs_data.K) // set q_l for allocated clusters   
                Weights_it(j,l) = std::log( gs_data.N(j,l) + gs_data.gamma[j] );
            
            else if(l==gs_data.K){ // set q_{K+1} for new cluster
                /*
                // OLD IMPLEMENTATION THAT IS WRONG
                Weights_it(j,l) =   log_const_new_cluster + 
                                    gs_data.d * std::log(gs_data.lambda) + 
                                    std::log(gs_data.gamma[j]) - 
                                    gs_data.gamma[j] * (gs_data.U[j] + 1.0);
                */
                Weights_it(j,l) =   log_const_new_cluster +
                                    std::log(gs_data.lambda * gs_data.gamma[j]) ;
                
            }
            else
                throw std::runtime_error("Error in GibbsSamplerMarginal::store_params_values, l can not be greater than gs_data.K ");
        }
    }
    out.log_q.push_back(Weights_it);
}



