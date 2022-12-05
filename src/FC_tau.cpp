#include "FC_tau.h"
#include <Rcpp.h>
#include <RcppEigen.h>


void FC_tau::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    //Retrive all data needed from gs_data
    const unsigned int& M = gs_data.M; // number of components
    const unsigned int& d = gs_data.d; // number of groups
    const unsigned int& K = gs_data.K; //number of clusters
    const std::vector<unsigned int>& n_j = gs_data.n_j; // number of observations per group
    const std::vector< std::vector<unsigned int>>& Ctilde = gs_data.Ctilde; // matrix of partition
    const std::vector<unsigned int>& N_k = gs_data.N_k;
    const std::vector<std::vector<double>>& data = gs_data.data; //matrix of data we don't copy it since data can be big but we use a pointer
    const std::string& prior = gs_data.prior; // identifier of the prior adopted for the model - togliamo la stringa e mettiamo una classe prior in modo che sia anche più leggibile
    // Initialize ind_i, ind_j
    std::vector<unsigned int> ind_i; // i index of C elements
    std::vector<unsigned int> ind_j;// j index of C elements

    
    //Rcpp::Rcout<<"Dentro FC_tau::Update"<<std::endl;
    //Rcpp::Rcout<<"gs_data.Mstar:"<<std::endl<<gs_data.Mstar<<std::endl;
    //Rcpp::Rcout<<"gs_data.M:"<<std::endl<<gs_data.M<<std::endl;
    //Rcpp::Rcout<<"gs_data.K:"<<std::endl<<gs_data.K<<std::endl;
    //Rcpp::Rcout<<"gs_data.mu.size():"<<std::endl<<gs_data.mu.size()<<std::endl;
    //Rcpp::Rcout<<"gs_data.sigma.size():"<<std::endl<<gs_data.sigma.size()<<std::endl;


    //Initialize tau according to new M
    gs_data.allocate_tau(gs_data.M); 

    if (prior == "Normal-InvGamma") {

        sample::rgamma Gamma;
        sample::rnorm rnorm;

        //NOT allocated tau

        for (unsigned int m = K; m < M; ++m){
             //Rcpp::Rcout<<"In marginal sampler case, the code should never reach this part"<<std::endl;
             double sigma2_na = 1 / Gamma(gs_engine, nu_0/2, 2/(nu_0 * sigma_0)); // Non allocated Components' variance
             double mu_na = rnorm(gs_engine, mu_0, std::sqrt(sigma2_na / k_0)); // Non allocated Components' mean
             gs_data.mu[m] = mu_na;
             gs_data.sigma[m] = sigma2_na;
             //Rcpp::Rcout << "Non Allocate: mu[" << m << "] = " << mu_na << std::endl;
             //Rcpp::Rcout << "sigma[" << m << "] = " << sigma2_na << std::endl;
        } 

        //Allocated tau
        for (unsigned int m = 0; m < K; ++m){
            ind_i.clear();
            ind_j.clear();
            for (unsigned int j = 0; j <d ; ++j) {
                for (unsigned int i = 0; i < n_j[j] ; ++i) {
                    //Rcpp::Rcout<<"Ctilde["<<j<<"]["<<i<<"]: "<<std::endl<<Ctilde[j][i]<<std::endl;
                    if(Ctilde[j][i] == m){
                        ind_i.push_back(i);
                        ind_j.push_back(j);
                    }
                }
            }
            //Esempio:
            /*
                
                Ctilde[0][0]: 0
                Ctilde[0][1]: 0
                Ctilde[1][0]: 0
                Ctilde[1][1]: 1

                Per m = 0 ho,
                Stampo ind_i: 0, 1, 0, 
                Stampo ind_j: 0, 0, 1, 

                Per m = 1 ho,
                Stampo ind_i: 1, 
                Stampo ind_j: 1, 

                Immagina Ctilde come una matrice, ind_i e ind_j mi dicono le coordinate che dovrò leggere componente per componente. Infatti nell'esempio,
                data[0,0] - data[0,1] - data[1,0] stanno nel primo cluster e data[1,1] sta nel secondo
            */

            //Rcpp::Rcout<<"Stampo ind_i: ";        
            //for(auto __v : ind_i)
                //Rcpp::Rcout<<__v<<", ";
            //Rcpp::Rcout<<std::endl;
//
            //Rcpp::Rcout<<"Stampo ind_j: ";        
            //for(auto __v : ind_j)
                //Rcpp::Rcout<<__v<<", ";
            //Rcpp::Rcout<<std::endl;

            //Rcpp::Rcout<<"N_k["<<m<<"]: "<<std::endl<<N_k[m]<<std::endl;
            double nu0_post = nu_0 + N_k[m];
            //Rcpp::Rcout<<"nu0_post:"<<std::endl<<nu0_post<<std::endl;
            double k0_post = k_0 + N_k[m];
            //Rcpp::Rcout<<"k0_post:"<<std::endl<<k0_post<<std::endl;
            double y_bar_clust= mean(ind_i,ind_j,data);
            //Rcpp::Rcout<<"y_bar_clust:"<<std::endl<<y_bar_clust<<std::endl;
            double s2_clust= var(y_bar_clust, ind_i, ind_j, data);
            //Rcpp::Rcout<<"s2_clust:"<<std::endl<<s2_clust<<std::endl;
            //if (is.na(s2_clust[k])){s2_clust[k] <- 0}
            double mu0_post = (k_0 * mu_0 + N_k[m] * y_bar_clust) / k0_post;
            //Rcpp::Rcout<<"mu0_post:"<<std::endl<<mu0_post<<std::endl;
            double sigma02_post =     1.0/(nu0_post) * ( nu_0*sigma_0 + 
                                                        ( (double)N_k[m] - 1.0 ) * s2_clust +
                                                        ( k_0 * (double)N_k[m] * (y_bar_clust - mu_0) * (y_bar_clust - mu_0) ) / (k0_post)
                                                       );
            //Rcpp::Rcout<<"sigma02_post:"<<std::endl<<sigma02_post<<std::endl;
            //Campionamento
            double sigma2_a = 1.0 / Gamma(gs_engine, nu0_post / 2.0, 2.0 / (nu0_post*sigma02_post) );
            double mu_a = rnorm(gs_engine, mu0_post, sqrt(sigma2_a / k0_post));
            gs_data.mu[m] = mu_a;
            gs_data.sigma[m] = sigma2_a;
            //Rcpp::Rcout << "Allocate: mu[" << m << "] = " << mu_a << std::endl;
            //Rcpp::Rcout << "sigma[" << m << "] = " << sigma2_a << std::endl;
        }
    }
}

// Function to compute the mean of the data (y_mean) for a group
double FC_tau::mean (const std::vector<unsigned int>& ind_i, const std::vector<unsigned int>& ind_j,
        const std::vector<std::vector<double>>& data){

    int count = 0;
    double sum = 0.0;
    for (size_t ii = 0; ii<ind_i.size(); ++ii) {
        sum += data.at(ind_j[ii]).at(ind_i[ii]);
        count++;
    }
    return sum/count;
}
// Function to compute the variance of the data
double FC_tau::var(double mean, const std::vector<unsigned int>& ind_i, const std::vector<unsigned int>& ind_j,
                    const std::vector<std::vector<double>>& data){
    if(ind_i.size() == 1)
        return 0.0;
    
    double vari=0.0;
    int count=0;
     for (size_t ii = 0; ii <ind_i.size() ; ++ii) {
         vari += (data.at(ind_j[ii]).at(ind_i[ii]) - mean) *(data.at(ind_j[ii]).at(ind_i[ii]) - mean) ;
         count++;
         //Rcpp::Rcout<<data.at(ind_j[ii]).at(ind_i[ii])<<",";
     }
     //Rcpp::Rcout<<std::endl;
     return vari/(count-1); 
}

