#include "FC_Partition_mv.h"
// Constructor 
Partition_mv::Partition_mv(std::string na, const unsigned int d, const std::vector<unsigned int> & n_j, bool FixPart) : FullConditional(na,FixPart), Partition_fixed(FixPart)
{
    //name = na;
    C.clear();
    for(size_t j = 0; j < d; j++){
        std::vector<unsigned int> row_j(n_j[j], 1);
        C.push_back(row_j);
    }
}

// update method
void Partition_mv::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    if(Partition_fixed){
        // Rcpp::Rcout << "Partition not updated because it is FIXED" << std::endl;
    }
    else{
        //Rcpp::Rcout<<" ----------  Dentro Partition  ----------"<<std::endl;
        // From gs_data all needed variable are retrived
        unsigned int k = gs_data.K; // number of cluster
        unsigned int d = gs_data.d; // number of group
        unsigned int M = gs_data.M; // number of components
        
        const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
        const std::vector<unsigned int>& n_j = gs_data.n_j;// number of observation per group
        const std::vector<double>& mu = gs_data.mu; // Vector of means
        const std::vector<double>& sigma = gs_data.sigma; // Vector of standard deviations
        const std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data;

        //Rcpp::Rcout<<"(k,Mstar,M) = "<<gs_data.K<<", "<<gs_data.Mstar<<", "<<gs_data.M<<std::endl;
//
        //Rcpp::Rcout<<"Stampo mu: ";        
        //for(auto __v : mu)
            //Rcpp::Rcout<<__v<<", ";
        //Rcpp::Rcout<<std::endl;
       //Rcpp::Rcout<<"Stampo sigma: ";        
        //for(auto __v : sigma)
            //Rcpp::Rcout<<__v<<", ";
        //Rcpp::Rcout<<std::endl;
        
 
        // Create vector to store probabilities for the M components
        GDFMM_Traits::VecRow probs_vec(M);
        // Initialization of probs_max
        double probs_max;
        //get sample index from GSL wrappers
        sample::sample_index sample_index;

        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){
                // compute "probability" each m component for y_ji 
                for(unsigned int m=0; m<M; m++){
                    //Rcpp::Rcout<<"S("<<j<<","<<m<<") = "<<S(j,m)<<std::endl;
                    //Rcpp::Rcout<<"log_dmvnorm(mv_data[j][i],mu[m],sigma[m]):"<<std::endl<<log_dmvnorm(mv_data[j][i],mu[m],sigma[m])<<std::endl;
                    //if(m>=K)
                        //Rcpp::Rcout<<"NON allocate: "<<std::endl;
                    
                    double likelihood_term{0.0};
                    if(gs_data.r > 0){
                        likelihood_term = log_dmvnorm2(mv_data[j][i],mu[m],mv_data[j][i].X_ji.transpose()*gs_data.beta.row(j), sigma[m]);
                    }
                    else
                        likelihood_term = log_dmvnorm2(mv_data[j][i],mu[m],sigma[m]);

                    if( std::abs(log_dmvnorm(mv_data[j][i],mu[m],sigma[m]) - likelihood_term) > 1e-6  ){
                        Rcpp::Rcout<<"i = "<<i<<std::endl;
                        Rcpp::Rcout<<"j = "<<j<<std::endl;
                        Rcpp::Rcout<<"mu[m]"<<mu[m]<<std::endl;
                        Rcpp::Rcout<<"sigma[m]"<<sigma[m]<<std::endl;
                        Rcpp::Rcout<<"log_dmvnorm(mv_data[j][i],mu[m],sigma[m]):"<<std::endl<<log_dmvnorm(mv_data[j][i],mu[m],sigma[m])<<std::endl;
                        Rcpp::Rcout<<"likelihood_term:"<<std::endl<<likelihood_term<<std::endl;
                        throw std::runtime_error("Error  in FC_Partition_mv, strano il termine di likelihood");
                    }
                    
                    probs_vec(m) = log(S(j,m)) + log_dmvnorm(mv_data[j][i],mu[m],sigma[m]);
                
                }
                // get the maximum "probability"
                probs_max = probs_vec.maxCoeff();

                // scale values of probs_vec
                for(unsigned int m=0; m<M; m++){
                    probs_vec(m) = exp(probs_vec(m) - probs_max);
                    //Rcpp::Rcout<<"probs_vec:"<<std::endl<<probs_vec<<std::endl;
                 //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                    if(std::isnan(probs_vec(m))){
                        Rcpp::Rcout<<"M = "<<M<<std::endl;
                        Rcpp::Rcout<<"nan m is m = "<<m<<std::endl;
                        Rcpp::Rcout<<"S("<<j<<","<<m<<") = "<<S(j,m)<<std::endl;
                        Rcpp::Rcout<<"mu["<<m<<"] = "<<mu[m]<<std::endl;
                        Rcpp::Rcout<<"sigma["<<m<<"] = "<<sigma[m]<<std::endl;
                        Rcpp::Rcout<<"log_dmvnorm(mv_data[j][i],mu[m],sigma[m]) = "<<log_dmvnorm(mv_data[j][i],mu[m],sigma[m])<<std::endl;
                        throw std::runtime_error("Error in Partition.cpp, get a nan in probs_vec ");
                    }
                }
                 //Rcpp::Rcout<<std::endl;
                // Assign y_ji to a component sampling from a multinomial
                C[j][i] = sample_index(gs_engine, probs_vec);
            }

        }

        // empty clust_out vector and set in order to reuse it
        clust_out.clear() ;
        s.clear();

        //Assign to each value of clust_out
        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){
                s.insert(C[j][i]); //insert every label inside a set
                clust_out.assign(s.begin(),s.end()); //get the vector of the label sorted and newly labeled e.g (0-1-2-3)
            }
        }

        k = clust_out.size(); //Set K=the size of clust out
        //Rcpp::Rcout<<"K = "<<k<<std::endl;
        gs_data.K = k; // updating K in the struct gs_data
        gs_data.allocate_N(k); // initialize N according to new K
        gs_data.update_Ctilde(C, clust_out);
                //Rcpp::Rcout<< "Numerosity in the "<< k << " clusters: ";
                //Rcpp::Rcout<<"Stampo gs_data.N_k: ";        
                //for(auto __v : gs_data.N_k)
                    //Rcpp::Rcout<<__v<<", ";
                //Rcpp::Rcout<<std::endl;

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


// density of N_{n_ji}(Y_ji; mu * ones(n_ji), var * diag(n_ji))
double Partition_mv::log_dmvnorm2(Individual data_ji, const double& mu, const double& var) const {
    double two_pi = 2.0 * M_PI;
    Eigen::Map<GDFMM_Traits::VecCol> eigen_data( &(data_ji.obs_ji[0]), data_ji.n_ji ); //cast observation into eigen form
    GDFMM_Traits::VecCol mu_term( GDFMM_Traits::VecCol::Constant(data_ji.n_ji, mu) ); // define a vector where each element is equal to mu
    GDFMM_Traits::VecCol mean( eigen_data - mu_term ); // compute the difference
    return(  -0.5*data_ji.n_ji*std::log(two_pi*var) - 
             (0.5/var) * mean.dot(mean) 
          );
}


// density of N_{n_ji}(Y_ji; mu * ones(n_ji) + X_ji*beta_j, var * diag(n_ji))
double Partition_mv::log_dmvnorm2(Individual data_ji, const double& mu,const GDFMM_Traits::VecCol& cov_term, const double& var) const {
    double two_pi = 2.0 * M_PI;
    Eigen::Map<GDFMM_Traits::VecCol> eigen_data( &(data_ji.obs_ji[0]), data_ji.n_ji ); //cast observation into eigen form
    GDFMM_Traits::VecCol mu_term( GDFMM_Traits::VecCol::Constant(data_ji.n_ji, mu) ); // define a vector where each element is equal to mu
    GDFMM_Traits::VecCol mean( eigen_data - mu_term - cov_term ); // compute the difference
    return(  -0.5*data_ji.n_ji*std::log(two_pi*var) - 
             (0.5/var) * mean.dot(mean) 
          );
}



// support method
double Partition_mv::log_dmvnorm(const Individual& data_ji, const double& mu, const double& var) const {
    double two_pi = 2.0 * M_PI;
    return(  -0.5*data_ji.n_ji*std::log(two_pi*var) - 
            (0.5/var) * ( (double)(data_ji.n_ji - 1)*data_ji.Vstar_ji + 
                          (double)data_ji.n_ji * (data_ji.Ybar_star_ji - mu)*(data_ji.Ybar_star_ji - mu) 
                        )  
          );
}
