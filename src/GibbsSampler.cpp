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

        // Manage cases if Partition is fixed
        std::vector<unsigned int> partition_vec;
        unsigned int Mstar0 = 0;

        if(Partition_fixed){
            partition_vec = Rcpp::as<std::vector<unsigned int>>(option["partition"]);
        }
        //else{
            Mstar0 = Rcpp::as<unsigned int>(option["Mstar0"]);
        //}
        // Read all hyper-parameters passed with option
        double nu = Rcpp::as<double>(option["nu"]);    
        double Lambda0 = Rcpp::as<double>(option["Lambda0"]);
        double mu0 = Rcpp::as<double>(option["mu0"]);
        double nu0 = Rcpp::as<double>(option["nu0"]);
        double sigma0 = Rcpp::as<double>(option["sigma0"]);
        double gamma0 = Rcpp::as<double>(option["gamma0"]);
        double h1 = Rcpp::as<double>(option["Adapt_MH_hyp1"]);
        double h2 = Rcpp::as<double>(option["Adapt_MH_hyp2"]);
        unsigned int pow = Rcpp::as<unsigned int>(option["Adapt_MH_power_lim"]);
        double adapt_var0 = Rcpp::as<double>(option["Adapt_MH_var0"]);
        double k0 = Rcpp::as<double>(option["k0"]);
        double a1 = Rcpp::as<double>(option["alpha_gamma"]);
        double b1 = Rcpp::as<double>(option["beta_gamma"]);
        double a2 = Rcpp::as<double>(option["alpha_lambda"]);
        double b2 = Rcpp::as<double>(option["beta_lambda"]);
        bool FixedU = !Rcpp::as<bool>(option["UpdateU"]);
        bool FixedM = !Rcpp::as<bool>(option["UpdateM"]);
        bool FixedGamma = !Rcpp::as<bool>(option["UpdateGamma"]);
        bool FixedS = !Rcpp::as<bool>(option["UpdateS"]);
        bool FixedTau = !Rcpp::as<bool>(option["UpdateTau"]);
        bool FixedLambda = !Rcpp::as<bool>(option["UpdateLambda"]);
        // Initialize gs_data with correct random seed, given Mstar and all data assigned to same cluster
        gs_data = GS_data( data, n_iter, burn_in, thin, random_engine,
                           Mstar0, Lambda0, mu0, nu0, sigma0, gamma0, P0_prior_name, partition_vec, nu);
        //Initialize Full Conditional Objects
        auto Partition_ptr = std::make_shared<Partition>("Partition", gs_data.d, gs_data.n_j, FixPart);
        auto Mstar_ptr = std::make_shared<FC_Mstar>("Mstar", FixPart, FixedM);
        auto gamma_ptr = std::make_shared<FC_gamma>("gamma", h1, h2, pow, adapt_var0, a1, b1, FixedGamma);
        auto tau_ptr = std::make_shared<FC_tau>("tau", nu0, sigma0, mu0, k0, FixedTau);
        auto U_ptr = std::make_shared<FC_U>("U", FixedU);
        auto S_ptr = std::make_shared<FC_S>("S", FixedS);
        auto lambda_ptr = std::make_shared<FC_Lambda>("lambda", a2, b2, FixedLambda);
        //Full Conditional vector that we will loop
        std::vector< std::shared_ptr<FullConditional> > fc{U_ptr,
                                                            Partition_ptr,
                                                            Mstar_ptr,
                                                            gamma_ptr,
                                                            S_ptr,
                                                            tau_ptr,
                                                            lambda_ptr
                                                            };
        //NOTE: it is enough to remove from fc the parametres that do not want to be updated. For example,
        /*
        // This removes Lambda and gamma
        std::vector< std::shared_ptr<FullConditional> > fc{U_ptr,
                                                            Partition_ptr,
                                                            Mstar_ptr,
                                                            S_ptr,
                                                            tau_ptr
                                                            };                                                                
        */
        std::swap(FullConditionals, fc);

        // Initialize return structure for S and tau
        out.S.reserve(burn_in + n_iter*thin);
        out.mu.reserve(burn_in + n_iter*thin);
        out.sigma.reserve(burn_in + n_iter*thin);
    }
}

void GibbsSampler::sample() {

    Progress progress_bar(burn_in + n_iter*thin, TRUE); // Initialize progress bar

    for(unsigned int it = 0; it <= burn_in + n_iter*thin; it++){
        // Sample from all full conditional
        this->GS_Step();

        //updating number of iterations necessary for MH algorithm
        gs_data.iterations = it;

        // If we are in the right iteration store needed values
        if(it>burn_in && (it-burn_in)%thin == 0){
            this->store_params_values();
        }

        if(!Partition_fixed){
            // Rcpp::Rcout<< "\nIn this iteration we obtain K: "<< gs_data.K << " M: " << gs_data.M
            //             <<"\n"<<std::endl;
        }

        progress_bar.increment(); //update progress bar
    }
}

void GibbsSampler::GS_Step() {
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
            gs_data.allocate_S(gs_data.M); // Initialize S according to new M 
            gs_data.allocate_tau(gs_data.M); // Initialize tau according to new M
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

void GibbsSampler::store_params_values() {
    
    out.K.push_back(gs_data.K); //modified
    out.Mstar.push_back(gs_data.Mstar); //modified
    if(!Partition_fixed){
        //out.K.push_back(gs_data.K);
        //out.Mstar.push_back(gs_data.Mstar);
        out.Ctilde.push_back(gs_data.Ctilde);
    }

    /*
    //POSSO TOGLIERE COMPLETAMENTE IL CALCOLO DEI w_jk??
    if(Partition_fixed){
        store_w_jk();
    }
    */

    // Common output values retrived
    //store_tau(); //Old version - ragazzi
    // Save tau
    out.mu.push_back( Rcpp::NumericVector (gs_data.mu.begin(),gs_data.mu.end()) );  //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
    out.sigma.push_back(  Rcpp::NumericVector (gs_data.sigma.begin(),gs_data.sigma.end())  ); //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
    // Save lambda - U - gamma
    out.lambda.push_back(gs_data.lambda);
    out.U.insert(out.U.end(), gs_data.U.begin(), gs_data.U.end() );
    out.gamma.insert(out.gamma.end(), gs_data.gamma.begin(), gs_data.gamma.end());

    //Update S
    out.S.push_back(gs_data.S);

    //Nuovo, salvo log_sum
    out.log_prod_psiU.push_back(gs_data.log_sum);
}

void GibbsSampler::store_w_jk(){
    
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

/* Old version - ragazzi
void GibbsSampler::store_tau(){
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
