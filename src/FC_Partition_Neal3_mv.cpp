#include "FC_Partition_Neal3_mv.h"

Partition_Neal3_mv::Partition_Neal3_mv( std::string _na, const unsigned int _d, const std::vector<unsigned int>& _n_j, bool _FixPart, 
                                         double _nu_0, double _sigma_0, double _mu_0, double _k_0 ):
                                         Partition_mv(_na,_d,_n_j,_FixPart), nu_0(_nu_0), sigma_0(_sigma_0), mu_0(_mu_0), k_0(_k_0){};

void Partition_Neal3_mv::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){

    if(Partition_fixed){
         //Rcpp::Rcout << "Partition not updated because it is FIXED" << std::endl;
    }
    else{
                //Rcpp::Rcout<<"*****************************"<<std::endl;
                //Rcpp::Rcout<<"Dentro Partition_Neal3_mv.cpp"<<std::endl;
        // From gs_data all needed variable are retrived
        const std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data; // get data structure and rename it for convenience
        unsigned int& K = gs_data.K; // number of cluster
        unsigned int& Mstar = gs_data.Mstar; // number of non active components
        const unsigned int& M = gs_data.M; // number of components
        const unsigned int& d = gs_data.d; // number of levels
        const unsigned int& r= gs_data.r; // number of regression coefficients
        const std::vector<unsigned int>& n_j = gs_data.n_j; // vector of length d, n_j[j] is the number of observed data in level j
        std::vector<unsigned int>& N_k = gs_data.N_k; // vector of length K, N_k[m] is the number of data assigned to cluster m
        std::vector<std::vector<unsigned int>>& Ctilde = gs_data.Ctilde; //vector of vector such that Ctilde[j][i] is the cluster membership for obs ji
        GDFMM_Traits::MatUnsCol& N = gs_data.N;   //dxK matrix, N(j,m) is the number of elements belonging to level j that are assigned to cluster m
        const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
        const GDFMM_Traits::MatRow& beta = gs_data.beta; // dxr matrix of regression coefficients
        GDFMM_Traits::vector_uset_uiui& cluster_indicies = gs_data.cluster_indicies; // vector with indicies for each cluster
        const bool& UseData = gs_data.UseData;
        
        // Quantities for computing marginal probabilities
        double nu0_post{0.0};
        double k0_post{0.0};
        double mu0_post{0.0};
        double sigma02_post{0.0};

        // Initialize ind_i, ind_j
        //std::vector<unsigned int> ind_i; // i index of C elements
        //std::vector<unsigned int> ind_j; // j index of C elements
        
        // Declare auxiliary quantities
        double log_probs_max;
        sample::sample_index sample_index;
        unsigned int new_Cji;
                
                /*
                Rcpp::Rcout<<"###################################"<<std::endl;
                Rcpp::Rcout<<"Stampo Ctilde INIZIALE: "<<std::endl;     
                for(auto __v : Ctilde){
                    for(auto ___v : __v){
                        Rcpp::Rcout<<___v<<", ";
                    }
                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<std::endl;
                Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                    Rcpp::Rcout<<__m<<": ";
                    std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                    [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                    Rcpp::Rcout<<std::endl;
                }
                */
                
        // Chinese Restaurant Franchise like process allocation
        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){

                // Remark: we are dealing with a CRF process, hence a table may be empty as long as at least one table serving the same
                // dish is occupied in one of the other levels/restaurants. Hence, N(j,m) may be 0 as long as N_k[m] > 0.
                // N_k[m] == 0 means that the global cluster m is now empty.

                    /*
                    Rcpp::Rcout<<"-----------------------------------"<<std::endl;
                    Rcpp::Rcout<<"(k,Mstar,M) = "<<gs_data.K<<", "<<gs_data.Mstar<<", "<<gs_data.M<<std::endl;
                    Rcpp::Rcout<<"Pre-update individuo "<<j<<", "<<i<<std::endl;
                    Rcpp::Rcout<<"Stampo N_k: ";        
                    for(auto __v : N_k)
                        Rcpp::Rcout<<__v<<", ";
                    Rcpp::Rcout<<std::endl;
    
                    Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;
                    
                    Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
                    for(auto __v : Ctilde){
                        for(auto ___v : __v){
                            Rcpp::Rcout<<___v<<", ";
                        }
                        Rcpp::Rcout<<std::endl;
                    }
                    Rcpp::Rcout<<std::endl;
                    
                    Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                    for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                        Rcpp::Rcout<<__m<<": ";
                        std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                        [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                        Rcpp::Rcout<<std::endl;
                    }
                    */
                    
    
                unsigned int& C_ji = Ctilde[j][i];  // what is the table assignment for obs ji? get its cluster membership
                        //Rcpp::Rcout<<"C_"<<j<<i<<":"<<std::endl<<C_ji<<std::endl;



                // remove obs ji from its cluster. In this step, both the local and the global counts must be updated as well as the sums in that cluster
                N_k[C_ji]--; // decrease the global counts
                N(j,C_ji)--; // decrease the local counts
                auto check_erase = cluster_indicies[C_ji].erase( std::make_pair(j,i) ); // erase element index from its cluster
                if(check_erase == 0){
                    Rcpp::Rcout<<"(j,i) = "<<j<<", "<<i<<std::endl;
                    Rcpp::Rcout<<"C_ji = "<<C_ji<<std::endl;
                    Rcpp::Rcout<<"check_erase = "<<check_erase<<std::endl;
                    Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                    for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                        Rcpp::Rcout<<__m<<": ";
                        std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                        [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                        Rcpp::Rcout<<std::endl;
                    }
                    Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;
                    throw std::runtime_error("Error, the element has not been eliminated from its cluster ");
                }

                // if the cluster becomes empty, then it must be removed. 
                // This is achived by replacing the now empty cluster with the last one.
                if(N_k[C_ji] == 0){
                    // last cluster to replace the now empty cluster. 
                    N_k[C_ji] = N_k[K-1];        // update the global counts
                    N_k.resize(K-1);             // eliminate the last cluster, now empty
                    N.col(C_ji) = N.col(K-1);    // update the local counts
                    cluster_indicies[C_ji] = cluster_indicies[K-1]; // update the vector of cluster indicies
                    cluster_indicies.resize(K-1); // eliminate the last cluster, now empty
                    N.conservativeResize(d,K-1); // eliminate the last cluster, now empty

                    // change all labels according to new labeling. 
                    for(unsigned int jj=0; jj<d; jj++){
                        for(unsigned int ii=0; ii<n_j[jj]; ii++){
                            if(Ctilde[jj][ii] == K-1)
                                Ctilde[jj][ii] = C_ji; // switch label for data currently in cluster K-1
                        }
                    }
                    
                    // In CRF I must set the current label for data ji to K-1, which is the now empty cluster. 
                    // This step is no longer necessary because I must compute the log probability of each possibile label for c_ji
                    //Ctilde[j][i] = K-1; 


                    K--; // decrese the current number of clusters 
                    Mstar++; // increase the number of non active components

                        /*
                        Rcpp::Rcout<<"***********************************"<<std::endl;
                        Rcpp::Rcout<<"Caso difficile - elimina il cluster"<<std::endl;
                        Rcpp::Rcout<<"Stampo N_k: ";        
                        for(auto __v : N_k)
                            Rcpp::Rcout<<__v<<", ";
                        Rcpp::Rcout<<std::endl;
    
                        Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;
    
                        Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
                        for(auto __v : Ctilde){
                            for(auto ___v : __v){
                                Rcpp::Rcout<<___v<<", ";
                            }
                            Rcpp::Rcout<<std::endl;
                        }
                        Rcpp::Rcout<<std::endl;
                        
                        Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                        for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                            Rcpp::Rcout<<__m<<": ";
                            std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                            [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                            Rcpp::Rcout<<std::endl;
                        }

                        Rcpp::Rcout<<"K = "<<K<<std::endl;
                        Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
                        */
                }


                GDFMM_Traits::VecRow log_probs_vec = GDFMM_Traits::VecRow::Constant(M, 0.0); // define vector of weights, currently in log-scale. lenght must be M

                // loop over current components 
                for(std::size_t m = 0; m < M; m++){
                            //Rcpp::Rcout<<"j = "<<j<<", i = "<<i<<", m = "<<m<<std::endl;
                    
                    // Place data ji in m-th component
                    Ctilde[j][i] = m; 

                    // m may exceede the size of cluster_indicies. If so, we are evaluating if (j,i) should form a cluster by its own
                    GDFMM_Traits::uset_uiui current_indicies; // Find data in component m

                    if(m < K){
                        current_indicies = cluster_indicies[m]; // copy all elements that are in such cluster
                    }
                    current_indicies.insert(std::make_pair(j,i)); // just add (j,i)

                    // Define useful quantities for posterior computations
                    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
                    double data_var_term{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
                    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
                    unsigned int sum_pi{0}; // sum(pi)
                            /*
                                // Find data in component m
                                // ---> This is no longer needed
                                ind_i.clear();
                                ind_j.clear();
                                for (unsigned int jj = 0; jj <d ; ++jj) {
                                    for (unsigned int ii = 0; ii < n_j[jj] ; ++ii) {
                                        //Rcpp::Rcout<<"Ctilde["<<j<<"]["<<i<<"]: "<<std::endl<<Ctilde[j][i]<<std::endl;
                                        if(Ctilde[jj][ii] == m){
                                            ind_i.push_back(ii);
                                            ind_j.push_back(jj);
                                        }
                                    }
                                }
                            */    

                            /*
                            Rcpp::Rcout<<"(K,Mstar,M) = "<<K<<", "<<Mstar<<", "<<M<<std::endl;
                            Rcpp::Rcout<<"m = "<<m<<std::endl;

                            Rcpp::Rcout<<"Stampo ind_j: ";        
                            for(auto __v : ind_j)
                                Rcpp::Rcout<<__v<<", ";
                            Rcpp::Rcout<<std::endl;
                            Rcpp::Rcout<<"Stampo ind_i: ";        
                            for(auto __v : ind_i)
                                Rcpp::Rcout<<__v<<", ";
                            Rcpp::Rcout<<std::endl;

                            Rcpp::Rcout<<"Stampo current_indicies "<<std::endl;
                            std::for_each(current_indicies.begin(), current_indicies.end(), 
                                            [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                            Rcpp::Rcout<<std::endl;
                            */
                            
                    //computed updated parameters
                    if(UseData){
                        if(r > 0)
                            std::tie(W,data_var_term,sum_piX2,sum_pi) = compute_cluster_summaries(current_indicies,mv_data,beta);
                        else
                            std::tie(W,data_var_term,sum_piX2,sum_pi) = compute_cluster_summaries(current_indicies,mv_data);
    
                                //Rcpp::Rcout<<"W:"<<std::endl<<W<<std::endl;
                                //Rcpp::Rcout<<"data_var_term:"<<std::endl<<data_var_term<<std::endl;
                                //Rcpp::Rcout<<"sum_piX2:"<<std::endl<<sum_piX2<<std::endl;
                                //Rcpp::Rcout<<"sum_pi:"<<std::endl<<sum_pi<<std::endl;
                    }

                    if(data_var_term < 0){
                        throw std::runtime_error("Error in Partition_Neal3_mv::update, data_var_term can not be negative ");
                    }
                    if( (sum_piX2 - (double)sum_pi*W*W) < -1e-8){
                        Rcpp::Rcout<<"W:"<<std::endl<<W<<std::endl;
                        Rcpp::Rcout<<"data_var_term:"<<std::endl<<data_var_term<<std::endl;
                        Rcpp::Rcout<<"sum_piX2:"<<std::endl<<sum_piX2<<std::endl;
                        Rcpp::Rcout<<"sum_pi:"<<std::endl<<sum_pi<<std::endl;
                        Rcpp::Rcout<<"sum_piX2 - (double)sum_pi*W*W:"<<std::endl<<sum_piX2 - (double)sum_pi*W*W<<std::endl;
                        Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                        for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                            Rcpp::Rcout<<__m<<": ";
                            std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                            [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                            Rcpp::Rcout<<std::endl;
                        }
                        throw std::runtime_error("Error in Partition_Neal3_mv::update, sum_piX2 - (double)sum_pi*W*W can not be negative ");
                    }

                    nu0_post = nu_0 + (double)sum_pi;
                    k0_post = k_0 + (double)sum_pi;
                    mu0_post = (k_0 * mu_0 + (double)sum_pi * W) / k0_post;
                    sigma02_post =     1.0/(nu0_post) * (   nu_0*sigma_0 +
                                                            data_var_term +
                                                            sum_piX2 - (double)sum_pi*W*W +
                                                            ( k_0*(double)sum_pi * (mu_0 - W) * (mu_0 - W) ) / (k0_post)
                                                        );

                            //Rcpp::Rcout<<"nu0_post:"<<std::endl<<nu0_post<<std::endl;
                            //Rcpp::Rcout<<"k0_post:"<<std::endl<<k0_post<<std::endl;
                            //Rcpp::Rcout<<"mu0_post:"<<std::endl<<mu0_post<<std::endl;
                            //Rcpp::Rcout<<"sigma02_post:"<<std::endl<<sigma02_post<<std::endl;
                    
                    // compute log weights
                    log_probs_vec[m] =  std::log( S(j,m) ); //prior term
                    if(UseData){
                        // add likelihood term
                        log_probs_vec[m] += log_I( mv_data[j][i].n_ji, mv_data[j][i].Ybar_star_ji, mv_data[j][i].Vstar_ji, mu0_post, k0_post, nu0_post, sigma02_post );
                    } 
                            //Rcpp::Rcout<<"std::log( S(j,m) ):"<<std::endl<<std::log( S(j,m) )<<std::endl;
                            //Rcpp::Rcout<<"log_I = "<<log_I( mv_data[j][i].n_ji, mv_data[j][i].Ybar_star_ji, mv_data[j][i].Vstar_ji,mu0_post, k0_post, nu0_post, sigma02_post )<<std::endl;
                            //Rcpp::Rcout<<"log_probs_vec[m]:"<<std::endl<<log_probs_vec[m]<<std::endl;
                    

                }
                
                        //Rcpp::Rcout<<"log_probs_vec:"<<std::endl<<log_probs_vec<<std::endl;
                
                // stable calculation of non-normalized weights in non-log scale
                log_probs_max = log_probs_vec.maxCoeff(); // get max value

                for(std::size_t m = 0; m < log_probs_vec.size(); m++){
                    log_probs_vec(m) = std::exp(log_probs_vec(m) - log_probs_max);
                     //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                    if(std::isnan(log_probs_vec(m)))
                        throw std::runtime_error("Error in Partition_Neal3_mv.cpp, get a nan in probs_vec ");
                }

                // Draw a sample of which cluster customer ji should belong to 
                new_Cji = sample_index(gs_engine, log_probs_vec); //values in log_probs_vec are no longer in log-scale here
                        //Rcpp::Rcout<<"------> new_C"<<j<<i<<":"<<std::endl<<new_Cji<<std::endl;


                // set a new cluster, if necessary
                //if( new_Cji == K) --> questo non ha piu senso, ora il cluster nuovo è qualunque new_Cji >= K. Devo riordinare in modo che siano sempre e solo i primi K valori
                // ad essere quelli rappresentativi delle componenti attive
                if( new_Cji >= K){
                            //Rcpp::Rcout<<"---> CREO UN NUOVO CLUSTER <---"<<std::endl;
                    N_k.resize(K+1);                         // allocate space in global counts for the new cluster
                    N_k[K] = 0;                              // values are assigned later
                    N.conservativeResize(d,K+1);             // allocate space in local counts for the new cluster
                    N.col(K) = VecUnsCol::Constant(d,0);     // values are assigned later. Note that here we are setting empty tables in all the other restaurants
                    cluster_indicies.resize(K+1);            // values are assigned later
                    new_Cji = K;                             // if here, new_Cji represents a non active component, hence it must be placed in the final cluster
                    K++;
                    Mstar--;
                }

                //Assign cluster membership and update counts
                Ctilde[j][i] = new_Cji; // set new label
                N_k[new_Cji]++;         // update global counts
                N(j,new_Cji)++;         // update local counts
                cluster_indicies[new_Cji].insert(std::make_pair(j,i)); // just add (j,i)
                    /*
                    Rcpp::Rcout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
                    Rcpp::Rcout<<"Post-update"<<std::endl;
                    Rcpp::Rcout<<"Stampo N_k: ";        
                    for(auto __v : N_k)
                        Rcpp::Rcout<<__v<<", ";
                    Rcpp::Rcout<<std::endl;
    
                    Rcpp::Rcout<<"N:"<<std::endl<<N<<std::endl;
    
                    Rcpp::Rcout<<"Stampo Ctilde: "<<std::endl;     
                    for(auto __v : Ctilde){
                        for(auto ___v : __v){
                            Rcpp::Rcout<<___v<<", ";
                        }
                        Rcpp::Rcout<<std::endl;
                    }
                    Rcpp::Rcout<<std::endl;

                    Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                    for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                        Rcpp::Rcout<<__m<<": ";
                        std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                        [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                        Rcpp::Rcout<<std::endl;
                    }
    
                    Rcpp::Rcout<<"K = "<<K<<std::endl;
                    Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
                    Rcpp::Rcout<<"+++++++++++++++++++++++++++++++++++"<<std::endl;
                    */
            }

        }

                /*
                Rcpp::Rcout<<"Stampo Ctilde FINALE: "<<std::endl;     
                for(auto __v : Ctilde){
                    for(auto ___v : __v){
                        Rcpp::Rcout<<___v<<", ";
                    }
                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"Stampo cluster_indicies "<<std::endl;  
                for (unsigned int __m = 0; __m < cluster_indicies.size(); ++__m){ // for each cluster
                    Rcpp::Rcout<<__m<<": ";
                    std::for_each(cluster_indicies[__m].cbegin(), cluster_indicies[__m].cend(), 
                                    [](const std::pair<unsigned int, unsigned int> p){Rcpp::Rcout<<"("<<p.first<<", "<<p.second<<") - ";});
                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<"###################################"<<std::endl;
                */

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

double Partition_Neal3_mv::log_I(unsigned int N_ji, double Ybar_star_ji, double Vstar_ji, double mu1, double k1, double nu1, double sigma1 ) const{
    double two_pi = 2.0 * M_PI;
    double nu1sigma1{nu1*sigma1};
    double coef1{ (k1 * (double)N_ji)/(k1 + (double)N_ji) };
    double term2{ ((double)N_ji - 1.0) * Vstar_ji};
    if(N_ji < 0)
        throw std::runtime_error("Error in log_I, N_ji can not be negative ");
    return{  0.5 * std::log(k1) - 
             0.5 * (double)N_ji * std::log(two_pi) + 
             0.5 * nu1 * std::log(0.5*nu1sigma1) + 
             std::lgamma(0.5*(nu1 + (double)N_ji)) - std::lgamma(0.5*nu1) - 
             0.5*(nu1 + (double)N_ji) * std::log( 0.5*( nu1sigma1 +  coef1*(mu1 - Ybar_star_ji)*(mu1 - Ybar_star_ji) + term2  ) )
         }; 
}


