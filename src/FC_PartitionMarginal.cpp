#include "FC_PartitionMarginal.h"

FC_PartitionMarginal::FC_PartitionMarginal( std::string _na, const unsigned int _d, const std::vector<unsigned int>& _n_j, bool _FixPart, 
                                            double _nu_0, double _sigma_0, double _mu_0, double _k_0 ):
                                            Partition(_na,_d,_n_j,_FixPart), nu_0(_nu_0), sigma_0(_sigma_0), mu_0(_mu_0), k_0(_k_0){};

// update method
void FC_PartitionMarginal::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    
    if(Partition_fixed){
         Rcpp::Rcout << "Partition not updated because it is FIXED" << std::endl;
    }
    else{
        //Rcpp::Rcout<<"Dentro Partition"<<std::endl;
        // From gs_data all needed variable are retrived
        const std::vector<std::vector<double>>& data = gs_data.data; //data
        const std::vector<std::vector<double>>& log_prob_marginal_data = gs_data.log_prob_marginal_data; //log marginal data

        //Rcpp::Rcout<<"Stampo log_prob_marginal_data: "<<std::endl;        
        //for(auto __v : log_prob_marginal_data){
            //for(auto ___v : __v){
                //Rcpp::Rcout<<___v<<", ";
            //}
            //Rcpp::Rcout<<std::endl;
        //}
        //Rcpp::Rcout<<std::endl;

        unsigned int& K = gs_data.K; // number of cluster
        const unsigned int& Lambda = gs_data.lambda; 
        const unsigned int& d = gs_data.d; // number of levels
        std::vector<unsigned int>& n_j = gs_data.n_j; // vector of length d, n_j[j] is the number of observed data in level j
        std::vector<unsigned int>& N_k = gs_data.N_k; // vector of length K, N_k[m] is the number of data assigned to cluster m
        std::vector< std::vector<unsigned int>>& Ctilde = gs_data.Ctilde; //vector of vector such that Ctilde[j][i] is the cluster membership for obs ji
        GDFMM_Traits::MatUnsCol& N = gs_data.N;   //dxK matrix, N(j,m) is the number of elements belonging to level j that are assigned to cluster m
        const std::vector<double>& U = gs_data.U; 
        const std::vector<double>& gamma = gs_data.gamma; 
        std::vector<double>& sum_cluster_elements = gs_data.sum_cluster_elements;
        std::vector<double>& squared_sum_cluster_elements = gs_data.squared_sum_cluster_elements;

        // Quantities for computing marginal probabilities
        const double& dof_prior{nu_0};
        const double& location_prior{mu_0};
        double dof_post{1.0};
        double location_post{0.0};
        double scale_post{1.0};

        // Declare auxiliary quantities
        double log_probs_max;
        sample::sample_index sample_index;
        unsigned int new_Cji;

        // Define the function to compute the constant that appears when computing the probability of creating a new cluster.
        // Note that, since U-gamma-Lambda are constants, this constant depends only on the number of clusters K. 
        // Hence its value must be updated only when K changes.
        auto compute_log_const_new_cluster = [&U, &gamma, &Lambda](const unsigned int& K){

            return (std::inner_product(U.cbegin(), U.cend(), gamma.cbegin(), 0.0, std::plus<>(),[&Lambda, &K](const double& U_j, const double& gamma_j)
                                                                                             {  
                                                                                                const double Psi_j{ 1/std::pow(1.0+U_j,gamma_j) }; // compute Psi_j
                                                                                                return (  -gamma_j*std::log(U_j+1.0) + 
                                                                                                          std::log((double)K+1.0+Lambda*Psi_j) - 
                                                                                                          std::log((double)K+Lambda*Psi_j)
                                                                                                       );
                                                                                             }
                                    )
                    );
        };
        double log_const_new_cluster{compute_log_const_new_cluster(K)};

        //Rcpp::Rcout<<"###################################"<<std::endl;
        // Chinese Restaurant Franchise process allocation
        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){

                        //Rcpp::Rcout<<"data["<<j<<"]["<<i<<"]:"<<std::endl<<data[j][i]<<std::endl;
                // Remark: we are dealing with a CRF process, hence a table may be empty as long as at least one table serving the same
                // dish is occupied in one of the other levels/restaurants. Hence, N(j,m) may be 0 as long as N_k[m] > 0.
                // N_k[m] == 0 means that the global cluster m is now empty.
                /*
                Rcpp::Rcout<<"-----------------------------------"<<std::endl;
                Rcpp::Rcout<<"Pre-update"<<std::endl;
                Rcpp::Rcout<<"Stampo N_k: ";        
                for(auto __v : N_k)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;

                Rcpp::Rcout<<"Stampo sum_cluster_elements: ";        
                for(auto __v : sum_cluster_elements)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"Stampo squared_sum_cluster_elements: ";        
                for(auto __v : squared_sum_cluster_elements)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
                for(auto __v : Ctilde){
                    for(auto ___v : __v){
                        Rcpp::Rcout<<___v<<", ";
                    }
                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"K = "<<K<<std::endl;

                Rcpp::Rcout<<"log_const_new_cluster:"<<std::endl<<log_const_new_cluster<<std::endl;
                */


                unsigned int C_ji = Ctilde[j][i];  // what is the table assignment for obs ji? get its cluster membership
                        //Rcpp::Rcout<<"C_"<<j<<i<<":"<<std::endl<<C_ji<<std::endl;

                // remove obs ji from its cluster. In this step, both the local and the global counts must be updated as well as the sums in that cluster
                N_k[C_ji]--; // decrease the global counts
                N(j,C_ji)--; // decrease the local counts
                sum_cluster_elements[C_ji] -= data[j][i];  //eliminate data_ji from the sum of elements in its cluster
                squared_sum_cluster_elements[C_ji] -= data[j][i]*data[j][i]; //eliminate data_ji from the sum of squared elements in its cluster

                // if the cluster becomes empty, then it must be removed. 
                // This is achived by replacing the now empty cluster with the last one.
                if(N_k[C_ji] == 0){
                    // last cluster to replace the now empty cluster. 
                    N_k[C_ji] = N_k[K-1];        // update the global counts
                    N_k.resize(K-1);             // eliminate the last cluster, now empty
                    N.col(C_ji) = N.col(K-1);    // update the local counts
                    N.conservativeResize(d,K-1); // eliminate the last cluster, now empty
                    sum_cluster_elements[C_ji] = sum_cluster_elements[K-1];                 // update the sum of variables
                    squared_sum_cluster_elements[C_ji] = squared_sum_cluster_elements[K-1]; // update the squared sum of variables
                    squared_sum_cluster_elements.resize(K-1);                               // eliminate the last cluster, now empty
                    sum_cluster_elements.resize(K-1);                                       // eliminate the last cluster, now empty

                    // change all labels according to new labeling. Data (there is only one and it is in position ji) with label C_ji is ruled out by setting its label equal to K-1
                    // while current data with label K-1 get label C_ji
                    
                    for(unsigned int jj=0; jj<d; jj++){
                        for(unsigned int ii=0; ii<n_j[jj]; ii++){
                            if(Ctilde[jj][ii] == K-1)
                                Ctilde[jj][ii] = C_ji; // switch label for data currently in cluster K-1
                        }
                    }
                    Ctilde[j][i] = K-1;

                    K--; // decrese the current number of clusters 
                    log_const_new_cluster = compute_log_const_new_cluster(K); // update constant

                    /*
                    Rcpp::Rcout<<"***********************************"<<std::endl;
                    Rcpp::Rcout<<"Caso difficile - elimina il cluster"<<std::endl;
                    Rcpp::Rcout<<"Stampo N_k: ";        
                    for(auto __v : N_k)
                        Rcpp::Rcout<<__v<<", ";
                    Rcpp::Rcout<<std::endl;

                    Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;

                    Rcpp::Rcout<<"Stampo sum_cluster_elements: ";        
                    for(auto __v : sum_cluster_elements)
                        Rcpp::Rcout<<__v<<", ";
                    Rcpp::Rcout<<std::endl;

                    Rcpp::Rcout<<"Stampo squared_sum_cluster_elements: ";        
                    for(auto __v : squared_sum_cluster_elements)
                        Rcpp::Rcout<<__v<<", ";
                    Rcpp::Rcout<<std::endl;

                    Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
                    for(auto __v : Ctilde){
                        for(auto ___v : __v){
                            Rcpp::Rcout<<___v<<", ";
                        }
                        Rcpp::Rcout<<std::endl;
                    }
                    Rcpp::Rcout<<std::endl;

                    Rcpp::Rcout<<"K = "<<K<<std::endl;

                    Rcpp::Rcout<<"log_const_new_cluster:"<<std::endl<<log_const_new_cluster<<std::endl;
                    */
                }


                GDFMM_Traits::VecRow log_probs_vec = GDFMM_Traits::VecRow::Constant(K+1, 0.0); // define vector of weights, currently in log-scale. lenght must be K+1

                // loop over current clusters 
                for(std::size_t l = 0; l < K; l++){
                    //computed updated parameters
                    dof_post = dof_prior + N_k[l];
                    double k_0_post = k_0 + N_k[l];
                    location_post = (k_0*location_prior + sum_cluster_elements[l])/k_0_post;
                    double sigma_0_post =   (1/dof_post) * (    (double)(N_k[l]-1) * gs_data.compute_var_in_cluster(l) +
                                                                dof_prior*sigma_0 +
                                                                (double)k_0/(k_0_post) * (double)N_k[l] * 
                                                                    (location_prior - sum_cluster_elements[l]/(double)N_k[l] ) * (location_prior - sum_cluster_elements[l]/(double)N_k[l] )
                                                            );

                    scale_post = std::sqrt(sigma_0_post*(k_0_post + 1.0)/k_0_post);
                            //Rcpp::Rcout<<"gs_data.compute_var_in_cluster("<<l<<") = "<<gs_data.compute_var_in_cluster(l)<<std::endl;
                            //Rcpp::Rcout<<"sigma_0_post = "<<sigma_0_post<<std::endl;
                            //Rcpp::Rcout<<"k_0_post = "<<k_0_post<<std::endl;
                            //Rcpp::Rcout<<"scale_post = "<<scale_post<<std::endl;
                    log_probs_vec[l] =  std::log( (double)N(j,l) + gamma[j] ) + 
                                        log_dnct(data[j][i], dof_post, location_post, scale_post);

                                        

                }
                
                // we are done looping over the already occupied clusters. Next, calculate the log probability of a new table.
                log_probs_vec[K] =  log_const_new_cluster + 
                                    std::log(gamma[j]) + d*std::log(Lambda) +
                                    log_prob_marginal_data[j][i];


                        //Rcpp::Rcout<<"log_probs_vec:"<<std::endl<<log_probs_vec<<std::endl;
                // stable calculation of non-normalized weights in non-log scale
                log_probs_max = log_probs_vec.maxCoeff(); // get max value

                for(std::size_t m = 0; m < log_probs_vec.size(); m++){
                    log_probs_vec(m) = std::exp(log_probs_vec(m) - log_probs_max);
                     //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                    if(std::isnan(log_probs_vec(m)))
                        throw std::runtime_error("Error in Partition.cpp, get a nan in probs_vec ");
                }

                // Draw a sample of which cluster customer ji should belong to 
                new_Cji = sample_index(gs_engine, log_probs_vec); //values in log_probs_vec are no longer in log-scale here
                        //Rcpp::Rcout<<"------> new_C"<<j<<i<<":"<<std::endl<<new_Cji<<std::endl;
                // set a new cluster, if necessary
                if( new_Cji == K){
                    N_k.resize(K+1);                         // allocate space in global counts for the new cluster
                    N_k[K] = 0;                              // values are assigned later
                    N.conservativeResize(d,K+1);             // allocate space in local counts for the new cluster
                    N.col(K) = VecUnsCol::Constant(d,0);     // values are assigned later. Note that here we are setting empty tables in all the other restaurants
                    sum_cluster_elements.resize(K+1);        // allocate space in vector having the sum of data in each cluster
                    squared_sum_cluster_elements.resize(K+1);// allocate space in vector having the squared sum of data in each cluster
                    K++;
                    log_const_new_cluster = compute_log_const_new_cluster(K); // update constant
                }

                //Assign cluster membership and update counts
                Ctilde[j][i] = new_Cji; // set new label
                N_k[new_Cji]++;         // update global counts
                N(j,new_Cji)++;         // update local counts
                sum_cluster_elements[new_Cji] += data[j][i];  //add data_ji to the sum of elements in its cluster
                squared_sum_cluster_elements[new_Cji] += data[j][i]*data[j][i]; //add data_ji to the sum of squared elements in its cluster

                /*
                Rcpp::Rcout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
                Rcpp::Rcout<<"Post-update"<<std::endl;
                Rcpp::Rcout<<"Stampo N_k: ";        
                for(auto __v : N_k)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;

                Rcpp::Rcout<<"Stampo sum_cluster_elements: ";        
                for(auto __v : sum_cluster_elements)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"Stampo squared_sum_cluster_elements: ";        
                for(auto __v : squared_sum_cluster_elements)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
                for(auto __v : Ctilde){
                    for(auto ___v : __v){
                        Rcpp::Rcout<<___v<<", ";
                    }
                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"K = "<<K<<std::endl;

                Rcpp::Rcout<<"log_const_new_cluster:"<<std::endl<<log_const_new_cluster<<std::endl;
                Rcpp::Rcout<<"+++++++++++++++++++++++++++++++++++"<<std::endl;
                */
            }

        }

        // In marginal sampler the non allocated components are not sampled. For sake of code generality, set M = K
        gs_data.M = K;

        if(K==0)
            throw std::runtime_error("K is 0, this should be impossible");
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


double FC_PartitionMarginal::dnct(const double& x, double const & n0, double const & mu0, const double& gamma0) const
{
    return( 1/gamma0 * gsl_ran_tdist_pdf( (x-mu0)/gamma0, n0 ) );
}

double FC_PartitionMarginal::log_dnct(const double& x, double const & n0, double const & mu0, const double& gamma0) const
{   
    if(n0 <= 0)
        throw std::runtime_error("Error in log_dnct, the degree of freedom must be strictly positive  ");
    if(gamma0 <= 0)
        throw std::runtime_error("Error in log_dnct, the scale must be strictly positive  ");
    return(   std::lgamma( (n0+1)/2 ) - std::lgamma( n0/2 ) - 0.5*std::log(M_PI*gamma0*n0) - 0.5*(n0 + 1)*std::log(1 + (1/n0)*(x-mu0)*(x-mu0)/(gamma0*gamma0)  )  );
}



// Stampa utile per dubugging
/*

Rcpp::Rcout<<"Stampo N_k: ";        
for(auto __v : N_k)
    Rcpp::Rcout<<__v<<", ";
Rcpp::Rcout<<std::endl;

Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;

Rcpp::Rcout<<"Stampo sum_cluster_elements: ";        
for(auto __v : sum_cluster_elements)
    Rcpp::Rcout<<__v<<", ";
Rcpp::Rcout<<std::endl;

Rcpp::Rcout<<"Stampo squared_sum_cluster_elements: ";        
for(auto __v : squared_sum_cluster_elements)
    Rcpp::Rcout<<__v<<", ";
Rcpp::Rcout<<std::endl;

Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
for(auto __v : Ctilde){
    for(auto ___v : __v){
        Rcpp::Rcout<<___v<<", ";
    }
    Rcpp::Rcout<<std::endl;
}
Rcpp::Rcout<<std::endl;

Rcpp::Rcout<<"K = "<<K<<std::endl;

Rcpp::Rcout<<"log_const_new_cluster:"<<std::endl<<log_const_new_cluster<<std::endl;
}

*/