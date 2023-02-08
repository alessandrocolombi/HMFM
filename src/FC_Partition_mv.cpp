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
        //Rcpp::Rcout<<"Dentro Partition"<<std::endl;
        // From gs_data all needed variable are retrived
        unsigned int k = gs_data.K; // number of cluster
        unsigned int d = gs_data.d; // number of group
        unsigned int M = gs_data.M; // number of components
        
        const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
        const std::vector<unsigned int>& n_j = gs_data.n_j;// number of observation per group
        const std::vector<double>& mu = gs_data.mu; // Vector of means
        const std::vector<double>& sigma = gs_data.sigma; // Vector of standard deviations
        const std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data;

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

                    //Rcpp::Rcout<<"log_dmvnorm(mv_data["<<j<<"]["<<i<<"],mu["<<m<<"],sigma["<<m<<"]) = "<<log_dmvnorm(mv_data[j][i],mu[m],sigma[m])<<std::endl;
                    probs_vec(m) = log(S(j,m)) + log_dmvnorm(mv_data[j][i],mu[m],sigma[m]);
                
                }
                // get the maximum "probability"
                probs_max = probs_vec.maxCoeff();

                // scale values of probs_vec
                for(unsigned int m=0; m<M; m++){
                    probs_vec(m) = exp(probs_vec(m) - probs_max);
                 //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                    if(std::isnan(probs_vec(m)))
                        throw std::runtime_error("Error in Partition.cpp, get a nan in probs_vec ");
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

// support method
double Partition_mv::log_dmvnorm(const Individual& data_ji, const double& mu, const double& var) const {
    double two_pi = 2.0 * M_PI;
    return(  -0.5*data_ji.n_ji*std::log(two_pi*var) - 
            (0.5/var) * ( (double)(data_ji.n_ji - 1)*data_ji.var_ji + 
                          (double)data_ji.n_ji * (data_ji.mean_ji - mu)*(data_ji.mean_ji - mu) 
                        )  
          );
}
