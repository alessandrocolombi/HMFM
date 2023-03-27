#include "Partition.h"
// Constructor
Partition::Partition(std::string na, const unsigned int d, const std::vector<unsigned int> & n_j,
                    bool FixPart) : FullConditional(na,FixPart), Partition_fixed(FixPart){
    //name = na;
    C.clear();
    for(size_t j = 0; j < d; j++){
        std::vector<unsigned int> row_j(n_j[j], 1);
        C.push_back(row_j);
    }
}

// update method
void Partition::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    if(Partition_fixed){
        // Rcpp::Rcout << "Partition not updated because it is FIXED" << std::endl;
    }
    else{
                Rcpp::Rcout<<"****************************"<<std::endl;
                Rcpp::Rcout<<"Dentro Partition"<<std::endl;
                Rcpp::Rcout<<"(M,K,Mstar) = "<<gs_data.M<<", "<<gs_data.K<<", "<<gs_data.Mstar<<std::endl;
        // From gs_data all needed variable are retrived
        unsigned int k = gs_data.K; // number of cluster
        unsigned int d = gs_data.d; // number of group
        unsigned int M = gs_data.M; // number of components

        //GDFMM_Traits::MatRow& S = gs_data.S; // MODIFIED, PUT IT CONST AGAIN
        const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
        const std::vector<unsigned int>& n_j = gs_data.n_j;// number of observation per group
        std::vector<double>& mu = gs_data.mu; // Vector of means
        std::vector<double>& sigma = gs_data.sigma; // Vector of standard deviations

                Rcpp::Rcout<<"Stampo mu: ";
                for(auto __v : mu)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

                Rcpp::Rcout<<"Stampo sigma: ";
                for(auto __v : sigma)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

        // Define data taken from gs_data
        const std::vector<std::vector<double>>& data = gs_data.data;
        // Create vector to store probabilities for the M components
        GDFMM_Traits::VecRow probs_vec(M);
        // Initialization of probs_max
        double probs_max;
        //get sample index from GSL wrappers
        sample::sample_index sample_index;

                Rcpp::Rcout<<"Prima di aggiornare, Ctilde:"<<std::endl;
                for(auto __v : gs_data.Ctilde){
                    for(auto __vv : __v)
                        Rcpp::Rcout<<__vv<<", ";

                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<std::endl;

        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){
                // compute "probability" each m component for y_ji
                for(unsigned int m=0; m<M; m++){

                    probs_vec(m) = log(S(j,m)) + log_norm(data[j][i], mu[m], sigma[m]);
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

                Rcpp::Rcout<<"Dopo aver aggiornato, C:"<<std::endl;
                for(auto __v : C){
                    for(auto __vv : __v)
                        Rcpp::Rcout<<__vv<<", ";

                    Rcpp::Rcout<<std::endl;
                }
                Rcpp::Rcout<<std::endl;

        // empty clust_out vector and set in order to reuse it
        clust_out.clear() ; // contiene i label originali. per esempio 1-4-5 se ho allocato osservazioni sono in corrispondenza dei valori unici che erano in posizione 1-4-5
        s.clear(); // s serve solo per creare clust_out

        //Assign to each value of clust_out
        for(unsigned int j=0; j<d; j++){
            for(unsigned int i=0; i<n_j[j]; i++){
                s.insert(C[j][i]); //insert every label inside a set
                clust_out.assign(s.begin(),s.end()); //get the vector of the label sorted and newly labeled e.g (0-1-2-3)
            }
        }

                Rcpp::Rcout<<"Stampo clust_out: ";
                for(auto __v : clust_out)
                    Rcpp::Rcout<<__v<<", ";
                Rcpp::Rcout<<std::endl;

        k = clust_out.size(); //Set K=the size of clust out
        //Rcpp::Rcout<<"K = "<<k<<std::endl;
        gs_data.K = k; // updating K in the struct gs_data
        gs_data.allocate_N(k); // initialize N according to new K
        gs_data.update_Ctilde(C, clust_out);

        /*
        I also need to update mu and sigma according to the new labels. Some of the old cluster may have disappeared,
        I must make sure that the first K components of mu and sigma represent the actual cluster parameters.
        NOTE that this operation may not be strictly necessary. Indeed, before updating the values of mu and sigma, the old
        values are erased. Before updating again the partition, such operation is always performed. Instead, this reordering
        is needed when such values are actually used, for example if some covariates are included in the model
        */
        std::vector<double> new_mu(gs_data.M, 0.0);
        std::vector<double> new_sigma(gs_data.M, 1.0);
        for(unsigned int m=0; m<k; m++){
            new_mu[m] = gs_data.mu[clust_out[m]];
            new_sigma[m] = gs_data.sigma[clust_out[m]];
        }
        // Now, set the non active components
        std::vector<unsigned int> old_numbering(gs_data.M);
        std::iota(old_numbering.begin(),old_numbering.end(),0);
        std::vector<unsigned int> diff;
        std::set_difference(old_numbering.begin(), old_numbering.end(), clust_out.begin(), clust_out.end(),
                                std::inserter(diff, diff.begin()));
                
        for(unsigned int m=k; m<gs_data.M; m++){
            new_mu[m] = gs_data.mu[diff[m-k]];
            new_sigma[m] = gs_data.sigma[diff[m-k]];
        }
        gs_data.mu = new_mu;
        gs_data.sigma = new_sigma;


                //Rcpp::Rcout<<"FINALE Dopo aver aggiornato, Ctilde:"<<std::endl;
                //for(auto __v : gs_data.Ctilde){
                    //for(auto __vv : __v)
                        //Rcpp::Rcout<<__vv<<", ";

                    //Rcpp::Rcout<<std::endl;
                //}
                //Rcpp::Rcout<<std::endl;

                //Rcpp::Rcout<<"Numerosita cluster, N_k"<<std::endl;
                //for(auto __v : gs_data.N_k)
                    //Rcpp::Rcout<<__v<<", ";
                //Rcpp::Rcout<<std::endl;
               //Rcpp::Rcout<<"Stampo mu: ";
                //for(auto __v : mu)
                    //Rcpp::Rcout<<__v<<", ";
                //Rcpp::Rcout<<std::endl;

                //Rcpp::Rcout<<"Stampo sigma: ";
                //for(auto __v : sigma)
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
double Partition::log_norm(double x, double u, double s) const {
    //const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;
    //return log((ONE_OVER_SQRT_2PI/std::sqrt(s))*std::exp(-0.5*(x-u)*(x-u)/s));
    double two_pi = 2.0 * M_PI;
    return(  -0.5*std::log(two_pi*s) -
            (0.5/s) * (x - u)*(x - u)
          );
}
