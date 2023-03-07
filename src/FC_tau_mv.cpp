#include "FC_tau_mv.h"
#include <Rcpp.h>
#include <RcppEigen.h>


void FC_tau_mv::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){

    //Rcpp::Rcout<<"Dentro FC_tau_mv.cpp "<<std::endl;
    //Retrive all data needed from gs_data
    const unsigned int& M = gs_data.M; // number of components
    const unsigned int& d = gs_data.d; // number of groups
    const unsigned int& r= gs_data.r; // number of regression coefficients
    const unsigned int& K = gs_data.K; //number of clusters
    const std::vector<unsigned int>& n_j = gs_data.n_j; // number of observations per group
    const std::vector< std::vector<unsigned int>>& Ctilde = gs_data.Ctilde; // matrix of partition
    const std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data; //matrix of data we don't copy it since data can be big but we use a pointer
    const std::string& prior = gs_data.prior; // identifier of the prior adopted for the model - togliamo la stringa e mettiamo una classe prior in modo che sia anche più leggibile
    const GDFMM_Traits::MatRow& beta = gs_data.beta; // dxr matrix of regression coefficients

    // Initialize ind_i, ind_j
    std::vector<unsigned int> ind_i; // i index of C elements
    std::vector<unsigned int> ind_j;// j index of C elements


    //Rcpp::Rcout<<"Dentro FC_tau_mv::Update"<<std::endl;
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
        double nu0_post{0.0};
        double k0_post{0.0};
        double mu0_post{0.0};
        double sigma02_post{0.0};
        double sigma2_m{0.0};
        double mu_m{0.0};

        double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
        double data_var_term{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
        double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
        unsigned int sum_pi{0}; // sum(pi)

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


            if(r > 0)
                std::tie(W,data_var_term,sum_piX2,sum_pi) = compute_cluster_summaries(ind_i,ind_j,mv_data,beta);
            else
                std::tie(W,data_var_term,sum_piX2,sum_pi) = compute_cluster_summaries(ind_i,ind_j,mv_data);

            //Rcpp::Rcout<<"W:"<<std::endl<<W<<std::endl;
            //Rcpp::Rcout<<"data_var_term:"<<std::endl<<data_var_term<<std::endl;
            //Rcpp::Rcout<<"sum_piX2:"<<std::endl<<sum_piX2<<std::endl;
            //Rcpp::Rcout<<"sum_pi:"<<std::endl<<sum_pi<<std::endl;

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

            //Campionamento
            sigma2_m = 1.0 / Gamma(gs_engine, nu0_post / 2.0, 2.0 / (nu0_post*sigma02_post) );
                    //sigma2_a = 1.0; // ONLY FOR DEBUGGING
                    //Rcpp::Rcout<<"Update component"<<std::endl;
                    //Rcpp::Rcout<<"Variance related coefficients"<<std::endl;
                    //Rcpp::Rcout<<"s2_clust = "<<std::endl<<s2_clust<<std::endl;
                    //Rcpp::Rcout<<"nu0_post = "<<nu0_post<<std::endl;
                    //Rcpp::Rcout<<"Prior term = "<<nu_0*sigma_0<<std::endl;
                    //Rcpp::Rcout<<"Variance term = "<<(double)N_k[m] - 1.0  * s2_clust<<std::endl;
                    //Rcpp::Rcout<<"squared difference term = "<<( k_0 * (double)N_k[m] * (y_bar_clust - mu_0) * (y_bar_clust - mu_0) ) / (k0_post)<<std::endl;
                    //Rcpp::Rcout<<"Expected value = "<< 
                                                        //nu_0*sigma_0/(nu0_post) <<" + "<<  ( (double)N_k[m] - 1.0 ) * s2_clust /nu0_post<< " + " << 
                                                        //( k_0 * (double)N_k[m] * (y_bar_clust - mu_0) * (y_bar_clust - mu_0) ) / (nu0_post*k0_post)<< " = " <<
                                                        //nu0_post*sigma02_post/nu0_post<<std::endl;
                    //Rcpp::Rcout<<"Mean related coefficients"<<std::endl;
                    //Rcpp::Rcout<<"n_k = "<<N_k[m]<<std::endl;
                    //Rcpp::Rcout<<"k_0 = "<<k_0<<std::endl;
                    //Rcpp::Rcout<<"peso 1: n_k/(n_k+k0) = "<< (double)N_k[m] / k0_post <<std::endl;
                    //Rcpp::Rcout<<"peso 2: k0/(n_k+k0) = "<< k_0 / k0_post <<std::endl;
                    //Rcpp::Rcout<<"y_bar_clust:"<<std::endl<<y_bar_clust<<std::endl;
                    //Rcpp::Rcout<<"mu_0:"<<std::endl<<mu_0<<std::endl;
                    //Rcpp::Rcout<<"mu0_post:"<<std::endl<<mu0_post<<std::endl;
            mu_m = rnorm(gs_engine, mu0_post, sqrt(sigma2_m / k0_post));
            gs_data.mu[m] = mu_m;
            gs_data.sigma[m] = sigma2_m;
            //Rcpp::Rcout << "Allocate: mu[" << m << "] = " << mu_m << std::endl;
            //Rcpp::Rcout << "sigma[" << m << "] = " << sigma2_m << std::endl;
        }

        //Rcpp::Rcout<<"Finito update tau"<<std::endl;
        //Rcpp::Rcout<<"Stampo gs_data.mu: ";        
        //for(auto __v : gs_data.mu)
            //Rcpp::Rcout<<__v<<", ";
        //Rcpp::Rcout<<std::endl;
//
        //Rcpp::Rcout<<"Stampo gs_data.sigma: ";        
        //for(auto __v : gs_data.sigma)
            //Rcpp::Rcout<<__v<<", ";
        //Rcpp::Rcout<<std::endl;
    }

}

std::tuple<double,double,double,unsigned int>  
FC_tau_mv::compute_cluster_summaries(  const std::vector<unsigned int>& ind_i, 
                                       const std::vector<unsigned int>& ind_j, 
                                       const std::vector<std::vector<Individual>>& data,
                                       const GDFMM_Traits::MatRow& beta )const
{
    // get number of elements in the cluster
    //unsigned int N_m{ind_i.size()};
    //Rcpp::Rcout<<"Dentro compute_cluster_summaries"<<std::endl;
    // initialize return objects
    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
    double mean_of_vars{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
    unsigned int sum_pi{0}; // sum(pi)

    // for each element in the cluster

    //  questa funzione non va bene se NON ho le covariate. il problema è l'ultimo termine, quello con beta e z.
    //  devo differenziare il calcolo di z in base a r=0 e r>0
    for(size_t ii = 0; ii <ind_i.size() ; ii++){

        const Individual& data_ji = data.at(ind_j[ii]).at(ind_i[ii]);
        
        //Rcpp::Rcout<<"data_ji.n_ji:"<<std::endl<<data_ji.n_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.Ybar_star_ji:"<<std::endl<<data_ji.Ybar_star_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.mean_ji:"<<std::endl<<data_ji.mean_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.z_ji:"<<std::endl<<data_ji.z_ji<<std::endl;

        sum_pi += data_ji.n_ji;
        mean_of_vars += (double)(data_ji.n_ji - 1)*data_ji.Ybar_star_ji;
        sum_piX2 += data_ji.n_ji * data_ji.Ybar_star_ji * data_ji.Ybar_star_ji;
        W += data_ji.n_ji*data_ji.mean_ji - data_ji.z_ji.dot( beta.row(ind_j[ii]) );

    }
    W = W/(double)sum_pi;

    //Rcpp::Rcout<<"W:"<<std::endl<<W<<std::endl;
    //Rcpp::Rcout<<"mean_of_vars:"<<std::endl<<mean_of_vars<<std::endl;
    //Rcpp::Rcout<<"sum_piX2:"<<std::endl<<sum_piX2<<std::endl;
    //Rcpp::Rcout<<"sum_pi:"<<std::endl<<sum_pi<<std::endl;
    
    return( std::tie(W,mean_of_vars,sum_piX2,sum_pi) );

}

std::tuple<double,double,double,unsigned int>  
FC_tau_mv::compute_cluster_summaries(  const std::vector<unsigned int>& ind_i, 
                                       const std::vector<unsigned int>& ind_j, 
                                       const std::vector<std::vector<Individual>>& data )const
{
    // get number of elements in the cluster
    //unsigned int N_m{ind_i.size()};
    //Rcpp::Rcout<<"Dentro compute_cluster_summaries"<<std::endl;
    // initialize return objects
    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
    double mean_of_vars{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
    unsigned int sum_pi{0}; // sum(pi)

    // for each element in the cluster

    //  questa funzione non va bene se NON ho le covariate. il problema è l'ultimo termine, quello con beta e z.
    //  devo differenziare il calcolo di z in base a r=0 e r>0
    for(size_t ii = 0; ii <ind_i.size() ; ii++){

        const Individual& data_ji = data.at(ind_j[ii]).at(ind_i[ii]);
        
        //Rcpp::Rcout<<"data_ji.n_ji:"<<std::endl<<data_ji.n_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.Ybar_star_ji:"<<std::endl<<data_ji.Ybar_star_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.mean_ji:"<<std::endl<<data_ji.mean_ji<<std::endl;
        sum_pi += data_ji.n_ji;
        mean_of_vars += (double)(data_ji.n_ji - 1)*data_ji.Ybar_star_ji;
        sum_piX2 += data_ji.n_ji * data_ji.Ybar_star_ji * data_ji.Ybar_star_ji;
        W += data_ji.n_ji*data_ji.mean_ji;

    }
    W = W/(double)sum_pi;

    //Rcpp::Rcout<<"W:"<<std::endl<<W<<std::endl;
    //Rcpp::Rcout<<"mean_of_vars:"<<std::endl<<mean_of_vars<<std::endl;
    //Rcpp::Rcout<<"sum_piX2:"<<std::endl<<sum_piX2<<std::endl;
    //Rcpp::Rcout<<"sum_pi:"<<std::endl<<sum_pi<<std::endl;
    
    return( std::tie(W,mean_of_vars,sum_piX2,sum_pi) );

}
/*
// Function to compute the mean of the data (y_mean) for a group
double FC_tau_mv::mean (const std::vector<unsigned int>& ind_i, const std::vector<unsigned int>& ind_j,
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
double FC_tau_mv::var(double mean, const std::vector<unsigned int>& ind_i, const std::vector<unsigned int>& ind_j,
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
*/
