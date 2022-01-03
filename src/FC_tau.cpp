#include "FC_tau.h"
#include "Partition.h"

void FC_tau::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    //Retrive all data needed from gs_data
    const unsigned int& Mna = gs_data.Mstar; // number of unallocated components
    const unsigned int& M = gs_data.M;
    const unsigned int& d = gs_data.d; // number of groups
    const unsigned int& K = gs_data.K; //number of clusters
    std::vector<unsigned int>& n_j= gs_data.n_j; // number of observations per group
    GDFMM_Traits::MatUnsCol& N = gs_data.N; // Matrix of observation per cluster per group
    std::vector< std::vector<unsigned int>>& Ctilde; // matrix of partition
    std::vector<unsigned int>& N_k = gs_data.N_k; //number of observations per cluster
    std::vector<std::vector<double>>& data=gs_data.data; //matrix of data we don't copy it since data can be big but we use a pointer
    Partition& p = gs_data.p; // actual partition of the data
    std::vector< std::vector<unsigned int>>& C=p.C; // C matrix of partition
    std::vector<unsigned int>& clust_out = p.clust_out; // Vector of clusters
    // RICONTROLLARE E CAPIRE DOVE METTERE CONST
    const char prior[]=gs_data.prior; // identifier of the prior adopted for the model
    // Initialize ind_i, ind_j
    std::vector<unsigned int> ind_i; // i index of C elements
    std::vector<unsigned int> ind_j;// j index of C elements

    if (std::strcmp(prior,"normal-inverse-gamma")) {

        sample::rgamma Gamma;
        sample::rnorm rnorm;

        //NOT allocated tau

        for (int m = K; m < M; ++m){
             double sigma2_na = 1 / Gamma(gs_engine, nu_0/2, (nu_0)*((sigma_0)/2)); // Non allocated Components' variance
             double mu_na = rnorm(gs_engine, mu_0, std::sqrt(sigma2_na / k_0)); // Non allocated Components' mean
             gs_data.mu[m] = mu_na;
             gs_data.sigma[m] = sigma2_na;
        }

        //Allocated tau

        for (int m = 0; m < K; ++m) {
            for (int j = 0; j <d ; ++j) {
                for (int i = 0; i < n_j[j] ; ++i) {
                    if(C[j][i] == clust_out[m]){
                        N(j,m)++;
                        ind_i.push_back(i);
                        ind_j.push_back(j);
                        Ctilde[j][m] = m;
                    }
                }
                N_k[m]+=N(j,m);
            }

            //set Ctilde in partition
            double nu_n_clust = nu_0 + N_k[m];
            double lpk = k_0 + N_k[m];
            double y_bar_clust= mean(ind_i, ind_j, data);
            double s2_clust= var(y_bar_clust, ind_i, ind_j, data);
            //if (is.na(s2_clust[k])){s2_clust[k] <- 0}
            double mu_n_clust = (k_0 * mu_0 + N_k[m] * y_bar_clust) / lpk;
            double sigma2_n_clust = (nu_0 * (sigma_0 * sigma_0) + (N_k[m] - 1) * s2_clust+ k_0 * N_k[m] * (y_bar_clust - mu_0) * (y_bar_clust - mu_0) / (lpk));

            //Campionamento
            double sigma2_a = 1 / Gamma(gs_engine, nu_n_clust/ 2, sigma2_n_clust / 2);
            double mu_a = rnorm(gs_engine, mu_n_clust, sqrt(sigma2_a / lpk));
            gs_data.mu[m] = mu_a;
            gs_data.sigma[m] = sigma2_a;
        }



}
}

// Function to compute the mean of the data (y_mean) for a group
double mean(std::vector<unsigned int> ind_i, std::vector<unsigned int> ind_j,
        const std::vector<std::vector<double>>& data){

    int count = 0;
    double sum = 0.0;
    for (int ii = 0; ii <ind_i.size(); ++ii) {
        sum += data.at(ind_j[ii]).at(ind_i[ii]);
        count++;
    }
    return sum/count;
}
// Function to compute the variance of the data
double var(double mean,std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j,
        const std::vector<std::vector<double>>& data){

    double var;
    int count=0;
     for (int ii = 0; ii <ind_i.size() ; ++ii) {
         var += (data.at(ind_j[ii]).at(ind_i[ii]) - mean) *
                (data.at(ind_j[ii]).at(ind_i[ii]) - mean) ;
         count++;
     }
     return var/count;
}

// mean e var vanno bene qui e poi si possono chiamare tranquillamente da dentro senza definire il namespace?
