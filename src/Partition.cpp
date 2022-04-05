#include "Partition.h"
// Constructor 
Partition::Partition(std::string na, const unsigned int d, const std::vector<unsigned int> & n_j,
                    bool FixPart) : Partition_fixed(FixPart){
    name = na;
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
        // From gs_data all needed variable are retrived
        unsigned int k = gs_data.K; // number of cluster
        unsigned int d = gs_data.d; // number of group
        unsigned int M = gs_data.M; // number of components
        
        //GDFMM_Traits::MatRow& S = gs_data.S; // MODIFIED, PUT IT CONST AGAIN
        const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
        const std::vector<unsigned int>& n_j = gs_data.n_j;// number of observation per group
        const std::vector<double>& mu = gs_data.mu; // Vector of means
        const std::vector<double>& sigma = gs_data.sigma; // Vector of standard deviations
        // Define data taken from gs_data
        const std::vector<std::vector<double>>& data = gs_data.data;
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

                    probs_vec(m) = log(S(j,m)) + log_norm(data[j][i], mu[m], sigma[m]);
                }
                // get the maximum "probability"
                probs_max = probs_vec.maxCoeff();

                // scale values of probs_vec
                for(unsigned int m=0; m<M; m++){
                    probs_vec(m) = exp(probs_vec(m) - probs_max);
                 //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                }
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
        gs_data.K = k; // updating K in the struct gs_data
        gs_data.allocate_N(k); // initialize N according to new K
        gs_data.update_Ctilde(C, clust_out);
        // Rcpp::Rcout<< "Numerosity in the "<< k << " clusters: ";
        /*
        for(unsigned int m=0; m<gs_data.K; m++){
            Rcpp::Rcout<<gs_data.N_k[m]<< " ";
        }
        */
    }
}

// support method
double Partition::log_norm(double x, double u, double s) const {
    const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;
    return log((ONE_OVER_SQRT_2PI/std::sqrt(s))*std::exp(-0.5*(x-u)*(x-u)/s));
}
