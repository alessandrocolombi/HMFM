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
        Rcpp::Rcout<<"Partition: "<<std::endl;
        for(unsigned int i = 0; i<partition_vec.size(); i++)
            Rcpp::Rcout<<partition_vec[i]<<"; ";
        
        Rcpp::Rcout<<std::endl;


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
        double k0 = Rcpp::as<double>(option["k0"]);
        double a1 = Rcpp::as<double>(option["alpha_gamma"]);
        double b1 = Rcpp::as<double>(option["beta_gamma"]);
        double a2 = Rcpp::as<double>(option["alpha_lambda"]);
        double b2 = Rcpp::as<double>(option["beta_lambda"]);
        bool FixedU = !Rcpp::as<bool>(option["UpdateU"]);
        bool FixedGamma = !Rcpp::as<bool>(option["UpdateGamma"]);
        bool FixedTau = !Rcpp::as<bool>(option["UpdateTau"]);
        bool FixedLambda = !Rcpp::as<bool>(option["UpdateLambda"]);

        // Initialize gs_data with correct random seed, given Mstar and all data assigned to same cluster
        gs_data = GS_data( data, n_iter, burn_in, thin, random_engine,
                           Mstar0, Lambda0, mu0, nu0, sigma0, gamma0, P0_prior_name, partition_vec);

        // Iinitialize sums in clusters
        gs_data.initialize_sums_in_clusters();
        gs_data.compute_log_prob_marginal_data(nu0, sigma0, mu0,  k0);

        //Initialize Full Conditional Objects
        auto PartitionMarginal_ptr = std::make_shared<FC_PartitionMarginal>("Partition", gs_data.d, gs_data.n_j, FixPart, nu0, sigma0, mu0, k0);
        auto gammaMarginal_ptr = std::make_shared<FC_gammaMarginal>("gamma", h1, h2, 10, gs_data.d, 1.0, a1, b1, FixedGamma);
        auto tau_ptr = std::make_shared<FC_tau>("tau", nu0, sigma0, mu0, k0, FixedTau);
        auto UMarginal_ptr = std::make_shared<FC_UMarginal>("U", FixedU, h1, h2, 10, gs_data.d, 1.0);
        auto lambdaMarginal_ptr = std::make_shared<FC_LambdaMarginal>("lambda", a2, b2, FixedLambda, h1, h2, 10, 1.0);

        //Full Conditional vector that we will loop
        std::vector< std::shared_ptr<FullConditional> > fc{ tau_ptr,
                                                            lambdaMarginal_ptr,
                                                            gammaMarginal_ptr,
                                                            UMarginal_ptr,
                                                            PartitionMarginal_ptr
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
        out.it_saved = 0;
    }
    else{
        throw std::runtime_error("Error, P0_prior_name must be equal to Normal-InvGamma. No other cases have been implemented.");
    }
}

void GibbsSamplerMarginal::sample() {

    Progress progress_bar(burn_in + n_iter*thin, TRUE); // Initialize progress bar
    for(unsigned int it = 0; it < burn_in + n_iter*thin; it++){

        // If we are in the right iteration store needed values
        // If burn_in is set to zero, the first values to be saved are the initial values.
        if(it >= burn_in && (it-burn_in)%thin == 0 && it < n_iter){
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

}



