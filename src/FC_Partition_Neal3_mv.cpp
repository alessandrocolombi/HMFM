#include "FC_Partition_Neal3_mv.h"

Partition_Neal3_mv::Partition_Neal3_mv( std::string _na, const unsigned int _d, const std::vector<unsigned int>& _n_j, bool _FixPart, 
                                         double _nu_0, double _sigma_0, double _mu_0, double _k_0 ):
                                         Partition_mv(_na,_d,_n_j,_FixPart), nu_0(_nu_0), sigma_0(_sigma_0), mu_0(_mu_0), k_0(_k_0){};

void Partition_Neal3_mv::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
/*
    if(Partition_fixed){
         //Rcpp::Rcout << "Partition not updated because it is FIXED" << std::endl;
    }
    else{
        Rcpp::Rcout<<"*****************************"<<std::endl;
        Rcpp::Rcout<<"Dentro Partition_Neal3_mv.cpp"<<std::endl;
        // From gs_data all needed variable are retrived
        const std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data; // get data structure and rename it for convenience
        unsigned int& K = gs_data.K; // number of cluster
        unsigned int& Mstar = gs_data.Mstar; // number of non active components
        const unsigned int& M = gs_data.M; // number of components
        const unsigned int& d = gs_data.d; // number of levels
        const std::vector<unsigned int>& n_j = gs_data.n_j; // vector of length d, n_j[j] is the number of observed data in level j
        std::vector<unsigned int>& N_k = gs_data.N_k; // vector of length K, N_k[m] is the number of data assigned to cluster m
        std::vector<std::vector<unsigned int>>& Ctilde = gs_data.Ctilde; //vector of vector such that Ctilde[j][i] is the cluster membership for obs ji
        GDFMM_Traits::MatUnsCol& N = gs_data.N;   //dxK matrix, N(j,m) is the number of elements belonging to level j that are assigned to cluster m
        const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights

        // Quantities for computing marginal probabilities
        double nu_1; 
        double sigma_1; 
        double mu_1; 
        double k_1;

        // Declare auxiliary quantities
        double log_probs_max;
        sample::sample_index sample_index;
        unsigned int new_Cji;

        //Rcpp::Rcout<<"###################################"<<std::endl;
        // Chinese Restaurant Franchise like process allocation
        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){

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

/*
                unsigned int C_ji = Ctilde[j][i];  // what is the table assignment for obs ji? get its cluster membership
                        //Rcpp::Rcout<<"C_"<<j<<i<<":"<<std::endl<<C_ji<<std::endl;

                // remove obs ji from its cluster. In this step, both the local and the global counts must be updated as well as the sums in that cluster
                N_k[C_ji]--; // decrease the global counts
                N(j,C_ji)--; // decrease the local counts
                //sum_cluster_elements[C_ji] -= data[j][i];  //eliminate data_ji from the sum of elements in its cluster
                //squared_sum_cluster_elements[C_ji] -= data[j][i]*data[j][i]; //eliminate data_ji from the sum of squared elements in its cluster

                // if the cluster becomes empty, then it must be removed. 
                // This is achived by replacing the now empty cluster with the last one.
                if(N_k[C_ji] == 0){
                    // last cluster to replace the now empty cluster. 
                    N_k[C_ji] = N_k[K-1];        // update the global counts
                    N_k.resize(K-1);             // eliminate the last cluster, now empty
                    N.col(C_ji) = N.col(K-1);    // update the local counts
                    N.conservativeResize(d,K-1); // eliminate the last cluster, now empty
                    //sum_cluster_elements[C_ji] = sum_cluster_elements[K-1];                 // update the sum of variables
                    //squared_sum_cluster_elements[C_ji] = squared_sum_cluster_elements[K-1]; // update the squared sum of variables
                    //squared_sum_cluster_elements.resize(K-1);                               // eliminate the last cluster, now empty
                    //sum_cluster_elements.resize(K-1);                                       // eliminate the last cluster, now empty

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
                    Mstar++;

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

/*
                }


                GDFMM_Traits::VecRow log_probs_vec = GDFMM_Traits::VecRow::Constant(M, 0.0); // define vector of weights, currently in log-scale. lenght must be M

                // loop over current components 
                for(std::size_t m = 0; m < M; l++){
                    //computed updated parameters
                    // SOME CODE HERE
                    log_probs_vec[m] =  std::log( S(j,m) ) + 
                                        log_I();

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
                if( new_Cji == K){
                    N_k.resize(K+1);                         // allocate space in global counts for the new cluster
                    N_k[K] = 0;                              // values are assigned later
                    N.conservativeResize(d,K+1);             // allocate space in local counts for the new cluster
                    N.col(K) = VecUnsCol::Constant(d,0);     // values are assigned later. Note that here we are setting empty tables in all the other restaurants
                    //sum_cluster_elements.resize(K+1);        // allocate space in vector having the sum of data in each cluster
                    //squared_sum_cluster_elements.resize(K+1);// allocate space in vector having the squared sum of data in each cluster
                    K++;
                    Mstar--;
                    //log_const_new_cluster = compute_log_const_new_cluster(K); // update constant
                }

                //Assign cluster membership and update counts
                Ctilde[j][i] = new_Cji; // set new label
                N_k[new_Cji]++;         // update global counts
                N(j,new_Cji)++;         // update local counts
                //sum_cluster_elements[new_Cji] += data[j][i];  //add data_ji to the sum of elements in its cluster
                //squared_sum_cluster_elements[new_Cji] += data[j][i]*data[j][i]; //add data_ji to the sum of squared elements in its cluster

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
/*

            }

        }


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

*/    
}