#include "ConditionalSampler.h"
#include <Rcpp.h>
#include <RcppEigen.h>

// Consrtructor
ConditionalSampler::ConditionalSampler( const Rcpp::List& _data_list, 
                                        unsigned int _n_it, unsigned int _b_in, unsigned int _thn,
                                        unsigned int _seed, std::string _P0_prior_name, bool _FixPart, 
                                        const Rcpp::List& _option) : random_engine(_seed), Partition_fixed(_FixPart), n_iter(_n_it),burn_in(_b_in),thin(_thn){
    // Extract hyper_parameter and initialization values from option
    if(_P0_prior_name == "Normal-InvGamma"){

        //Rcpp::Rcout<<"ConditionalSampler.cpp - Initialization"<<std::endl;

        // Manage cases if Partition is fixed
        std::vector<unsigned int> partition_vec{ Rcpp::as<std::vector<unsigned int>>(_option["partition"]) }; 
        unsigned int Mstar0{Rcpp::as<unsigned int>(_option["Mstar0"])}; 
        
        // Stampe
        //Rcpp::Rcout<<"Mstar0 = "<<Mstar0<<std::endl;
        //Rcpp::Rcout<<"Partition: "<<std::endl;
        //for(unsigned int i = 0; i<partition_vec.size(); i++)
            //Rcpp::Rcout<<partition_vec[i]<<"; ";
//        
        //Rcpp::Rcout<<std::endl;
        //partition_vec = Rcpp::as<std::vector<unsigned int>>(_option["partition"]);
        // Read data passed with _data_list
        n = Rcpp::as<unsigned int>(_data_list["n"]);
        d = Rcpp::as<unsigned int>(_data_list["d"]);
        r = Rcpp::as<unsigned int>(_data_list["r"]);
        n_j =  Rcpp::as<std::vector<unsigned int>>(_data_list["n_j"]);
        ID_i = Rcpp::as<std::vector<std::string>>(_data_list["ID_i"]);
        //s_i =  Rcpp::as<std::vector<unsigned int>>(_data_list["s_i"]);
        
        Rcpp::IntegerMatrix N_ji    = Rcpp::as<Rcpp::IntegerMatrix>(_data_list["N_ji"]);
        Rcpp::NumericMatrix mean_ji = Rcpp::as<Rcpp::NumericMatrix>(_data_list["mean_ji"]);
        Rcpp::NumericMatrix var_ji  = Rcpp::as<Rcpp::NumericMatrix>(_data_list["var_ji"]);
        Rcpp::List obs = Rcpp::as<Rcpp::List>(_data_list["observations"]);
        Rcpp::List covariates = Rcpp::as<Rcpp::List>(_data_list["covariates"]);
        bool IncludeCovariates = Rcpp::as<bool>(_option["IncludeCovariates"]);
        std::vector<std::vector<Individual>> data;
        data.resize(d);
        for(size_t j = 0; j < d; j++)
            data[j].reserve(n_j[j]);

        //Rcpp::Rcout<<"data.size() = "<<data.size()<<std::endl;
        for(size_t j = 0; j < d; j++){
            Rcpp::List obs_j = obs[j]; // get list of all individualus at level j
            Rcpp::List X_j;   
            if(IncludeCovariates)
                X_j = covariates[j]; // get list of all individualus at level j
            for(size_t i = 0; i < n; i++){
                if(N_ji(j,i) > 0){ // add only individuals who have at least one observation
                            //Rcpp::Rcout<<"-----------------------------------"<<std::endl;
                            //Rcpp::Rcout<<"Creo livello "<<j<<", individuo "<<i<<std::endl;
                    std::vector<double> obs_ji = Rcpp::as<std::vector<double>>(obs_j[i]); // get data for individual i in level j
                    if(IncludeCovariates){
                        Rcpp::NumericMatrix X_ji = Rcpp::as<Rcpp::NumericMatrix>(X_j[i]); // get design matrix for individual i in level j
                        Individual data_ji( ID_i[i], (unsigned int)N_ji(j,i), mean_ji(j,i), var_ji(j,i), obs_ji, X_ji );
                        data[j].push_back(data_ji); // add individual to the list
                    }
                    else{
                        Individual data_ji( ID_i[i], (unsigned int)N_ji(j,i), mean_ji(j,i), var_ji(j,i), obs_ji ); // get design matrix for individual i in level j
                        data[j].push_back(data_ji);   // add individual to the list
                    }
                }

            }
            //Rcpp::Rcout<<"data["<<j<<"].size() = "<<data[j].size()<<std::endl;
        }

                //Rcpp::Rcout<<"data.size() = "<<data.size()<<std::endl;
                //for(size_t j = 0; j < d; j++)
                    //for(size_t i = 0; i < data[j].size(); i++)
                        //Rcpp::Rcout<<"data["<<j<<","<<i<<"] = "<<data[j][i].mean_ji<<std::endl;


        // Read all hyper-parameters passed with _option
        double Lambda0 = Rcpp::as<double>(_option["Lambda0"]);
        double mu0 = Rcpp::as<double>(_option["mu0"]);
        double nu0 = Rcpp::as<double>(_option["nu0"]);
        double sigma0 = Rcpp::as<double>(_option["sigma0"]);
        double gamma0 = Rcpp::as<double>(_option["gamma0"]);

        // Read Rcpp objects and transform them in Eigen objects
        Eigen::Map<GDFMM_Traits::MatCol> Sigma0_temp = Rcpp::as<Eigen::Map<GDFMM_Traits::MatCol>>( Rcpp::as<Rcpp::NumericMatrix>(_option["Sigma0"]) ); 
        Eigen::Map<GDFMM_Traits::VecCol> beta0_temp  = Rcpp::as<Eigen::Map<GDFMM_Traits::VecCol>>( Rcpp::as<Rcpp::NumericVector>(_option["beta0"]) ); 
        GDFMM_Traits::VecCol beta0 = beta0_temp; 
        GDFMM_Traits::MatCol Sigma0 = Sigma0_temp; 

        double h1 = Rcpp::as<double>(_option["Adapt_MH_hyp1"]);
        double h2 = Rcpp::as<double>(_option["Adapt_MH_hyp2"]);
        unsigned int pow = Rcpp::as<unsigned int>(_option["Adapt_MH_power_lim"]);
        unsigned int proposal_Mstar = Rcpp::as<unsigned int>(_option["proposal_Mstar"]);
        double adapt_var0 = Rcpp::as<double>(_option["Adapt_MH_var0"]);
        double k0 = Rcpp::as<double>(_option["k0"]);
        double a1 = Rcpp::as<double>(_option["alpha_gamma"]);
        double b1 = Rcpp::as<double>(_option["beta_gamma"]);
        double a2 = Rcpp::as<double>(_option["alpha_lambda"]);
        double b2 = Rcpp::as<double>(_option["beta_lambda"]);
        std::vector<double> init_mean_clus{ Rcpp::as<std::vector<double>>(_option["init_mean_cluster"]) };
        std::vector<double> init_var_clus{ Rcpp::as<std::vector<double>>(_option["init_var_cluster"]) };
        bool FixedU = !Rcpp::as<bool>(_option["UpdateU"]);
        bool FixedM = !Rcpp::as<bool>(_option["UpdateM"]);
        bool FixedGamma = !Rcpp::as<bool>(_option["UpdateGamma"]);
        bool FixedS = !Rcpp::as<bool>(_option["UpdateS"]);
        bool FixedTau = !Rcpp::as<bool>(_option["UpdateTau"]);
        bool FixedLambda = !Rcpp::as<bool>(_option["UpdateLambda"]);
        bool FixedBeta = !Rcpp::as<bool>(_option["UpdateBeta"]);
        // Initialize gs_data with correct random seed, given Mstar and all data assigned to same cluster
        gs_data = GS_data(  data, n_j, d, r, random_engine, 
                            Mstar0, Lambda0, mu0, nu0, sigma0, gamma0, 
                            beta0, Sigma0, 
                            init_mean_clus, init_var_clus, _P0_prior_name, 
                            partition_vec  );

        //Initialize Full Conditional Objects
        auto Partition_ptr = std::make_shared<Partition_mv>("Partition", d, n_j, _FixPart);
        auto Mstar_ptr = std::make_shared<FC_Mstar>("Mstar", proposal_Mstar, _FixPart, FixedM);
        auto gamma_ptr = std::make_shared<FC_gamma>("gamma", h1, h2, pow, d, adapt_var0, a1, b1, FixedGamma);
        auto tau_ptr = std::make_shared<FC_tau_mv>("tau", nu0, sigma0, mu0, k0, FixedTau);
        auto U_ptr = std::make_shared<FC_U>("U", FixedU);
        auto S_ptr = std::make_shared<FC_S>("S", FixedS);
        auto lambda_ptr = std::make_shared<FC_Lambda>("lambda", a2, b2, FixedLambda);
        auto beta_ptr = std::make_shared<FC_beta_mv>("beta", beta0, Sigma0, FixedBeta);


        //Full Conditional vector that we will loop
        std::vector< std::shared_ptr<FullConditional> > fc_nocov{tau_ptr,
                                                            S_ptr,
                                                            lambda_ptr,
                                                            Partition_ptr,
                                                            U_ptr,
                                                            gamma_ptr,
                                                            Mstar_ptr
                                                            };

        //Full Conditional vector that we will loop
        std::vector< std::shared_ptr<FullConditional> > fc_cov{beta_ptr,
                                                                tau_ptr,
                                                                S_ptr,
                                                                lambda_ptr,
                                                                Partition_ptr,
                                                                U_ptr,
                                                                gamma_ptr,
                                                                Mstar_ptr
                                                            }; 

        if(IncludeCovariates)            
            std::swap(FullConditionals, fc_cov);
        else
            std::swap(FullConditionals, fc_nocov);    

        // Initialize return structure 
        out.S.reserve(n_iter);
        out.mu.reserve(n_iter);
        out.sigma.reserve(n_iter);
        
        // Initialize beta return structure - only if needed
        if(IncludeCovariates)  
            out.beta.reserve(n_iter);

    }
    else{
        throw std::runtime_error("Error, P0_prior_name must be equal to Normal-InvGamma. No other cases have been implemented.");
    }
}

void ConditionalSampler::sample() {

    Progress progress_bar(burn_in + n_iter*thin, TRUE); // Initialize progress bar
    for(unsigned int it = 0; it < burn_in + n_iter*thin; it++){

        // If we are in the right iteration store needed values
        // If burn_in is set to zero, the first values to be saved are the initial values.
        if(it>=burn_in && (it-burn_in)%thin == 0){
            //Rcpp::Rcout<<"it = "<<it<<std::endl;
            this->store_params_values();
        }

        // Sample from all full conditional
        this->GS_Step();

        //updating number of iterations necessary for MH algorithm
        gs_data.iterations = it;

        progress_bar.increment(); //update progress bar
    }
    
}

void ConditionalSampler::GS_Step() {
    
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
            //Rcpp::Rcout<<"Aggiorno M senza aggiornare Mstar"<<std::endl;
            gs_data.M = gs_data.K + gs_data.Mstar;
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

void ConditionalSampler::store_params_values() {
    
    out.K.push_back(gs_data.K); //modified
    out.Mstar.push_back(gs_data.Mstar); //modified
    if(!Partition_fixed){
        //out.K.push_back(gs_data.K);
        //out.Mstar.push_back(gs_data.Mstar);
        out.Ctilde.push_back(gs_data.Ctilde);
    }

    // Common output values retrived
    //store_tau(); //Old version - ragazzi
    // Save tau
    out.mu.push_back( Rcpp::NumericVector (gs_data.mu.begin(),gs_data.mu.end()) );  //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
    out.sigma.push_back(  Rcpp::NumericVector (gs_data.sigma.begin(),gs_data.sigma.end())  ); //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
    // Save lambda - U - gamma
    out.lambda.push_back(gs_data.lambda);
    out.U.insert(out.U.end(), gs_data.U.begin(), gs_data.U.end() );
    out.gamma.insert(out.gamma.end(), gs_data.gamma.begin(), gs_data.gamma.end());

    //Save S
    out.S.push_back(gs_data.S);

    //Save beta - if needed
    if(gs_data.r > 0)  
        out.beta.push_back(gs_data.beta);

    //Save log sum to perform debugging/posterior checks
    out.log_prod_psiU.push_back(gs_data.log_sum);
    
}



/* Old version - ragazzi

void ConditionalSampler::store_w_jk(){
    
    unsigned int current_it = (gs_data.iterations - burn_in)/thin;
    unsigned int d = gs_data.S.rows();
    unsigned int K = gs_data.S.cols();
    
    if(out.w_jk.empty()){
        GDFMM_Traits::MatRow w_j(K, n_iter);
        for(unsigned int j = 0; j < d; j++){
            out.w_jk.push_back(w_j);
        }
    }

    Eigen::VectorXd T = gs_data.S.rowwise().sum();

    for(unsigned int j = 0; j < d; j++){
        for(unsigned int k = 0; k < K; k++){
            out.w_jk[j](k, current_it - 1) = gs_data.S(j,k)/T(j);
        }
    }
}

void ConditionalSampler::store_tau(){
    unsigned int current_it = (gs_data.iterations - burn_in)/thin;
    unsigned int current_K = gs_data.K;
    unsigned int size_tau = out.mu.size();

    // Check that current_K doesn't exceed size_tau
    if(out.mu.empty()){
        for(unsigned int k = 0; k < current_K; k++){
            out.mu.push_back(std::vector<double>{gs_data.mu[k]});  
            out.sigma.push_back(std::vector<double>{gs_data.sigma[k]});
        }
        return;
    }
    
    if(current_K > size_tau){
        // resize output tau accordingly and initialize values for past iterations
        out.mu.resize(current_K);
        out.sigma.resize(current_K);
        for(unsigned int i = size_tau; i < current_K; i++){
            out.mu[i] = std::vector<double>(current_it-1, std::nan("") );
            out.sigma[i] = std::vector<double>(current_it-1, std::nan("") );
        }
    }

    // Update tau with current value
    for(unsigned int k = 0; k < current_K; k++){
        out.mu[k].push_back(gs_data.mu[k]);
        out.sigma[k].push_back(gs_data.sigma[k]);
    }

    // Fill other tau, if there exist
    if(current_K < size_tau){
        for(unsigned int i = current_K; i < size_tau; i++){
            out.mu[i].push_back(std::nan("") );
            out.sigma[i].push_back(std::nan("") );
        }
    }
}
*/
