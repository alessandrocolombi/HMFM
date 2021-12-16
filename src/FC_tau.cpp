
#include "FC_tau.h"



void FC_tau::update(GS_data& gs_data, sample::GSL_RNG gs_engine, const string &c, Partition& p) {
    unsigned int Mna = gs_data.Mstar;
    unsigned int d = gs_data.d;
    unsigned int K=gs_data.K;
    std::vector<unsigned int> n_j= gs_data.n_j;
    std::vector< std::vector<unsigned int>> C=p.C;
    std::vector<unsigned int> clust_out= p.clust_out;
    GDFMM_Traits::MatUnsCol N=gs_data.N;
    std::vector<unsigned int> ind_i;
    std::vector<unsigned int> ind_j;
    std::vector< std::vector<unsigned int>> Ctilde=p.Ctilde;

    if (c == 'normal-inverse-gamma') {
        //tau-nonallocate
        sample::rgamma Gamma;
        sample::rnorm rnorm;
        for (int m=0;m<Mna;++m){
             float sigma2_na=1 / Gamma(gs_engine, nu_0/2, (nu_0)*((sigma_0)/2)); // vardellecomponenti
             float mu_na=rnorm( mu_0, std::sqrt(sigma2_na / k_0)); // mediadellecomponenti

             GS_data::mu.pushback(mu_na);
             GS_data::sigma.pushback(sigma2_na);
        }
         //tau-allocate
        for (int k= 0; k <K; ++k) {
            for (int j = 0; j <d ; ++j) {
                for (int i = 0; i <n_j[j] ; ++i) {
                    if(C[j][i]==clust_out[k]){
                        N(j,k)++;
                        ind_i.pushback(i);
                        ind_j.pushback(j);
                        Ctilde[j][k]=k;
                    }
                }
            }

            nu_n_clust = nu_0 + N_k[k];

            lpk = k_0 + N_k[k];
            y_bar_clust= mean(ind_i,ind_j);
            s2_clust= var(mean,ind_i, ind_j);
//if (is.na(s2_clust[k])){s2_clust[k] <- 0}
            mu_n_clust = (k_0 * mu_0 + N_k[k] * y_bar_clust) / lpk;
            sigma2_n_clust = (nu_0 * (sigma_0 ^ 2) + (N_k[k] - 1) * s2_clust+ k_0 * N_k[k]  *(y_bar_clust - mu_0) ^ 2 / (lpk));

            //Campionamento
            sigma2_a <- 1 / rgamma(gs_engine, nu_n_clust/ 2, sigma2_n_clust / 2);
            mu_a<- rnorm(gs_engine, mu_n_clust, sqrt(sigma2_a / lpk));
            tau_a[k,] <- c(mu_a[k], sigma2_a[k]);
        }



}
}

double mean(std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j){
    int count=0;
    double sum=0;
    for (int m = 0; m <ind_i ; ++m) {
        sum+=gs_data:data[ind_j[m]][ind_i[m]];
        count++;

    }
    return sum/count;
}
 double var(double mean,std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j){
    double var
    int count=0
     for (int m = 0; m <ind_i ; ++m) {
         var+=(gs_data:data[ind_j[m]][ind_i[m]] - mean) * (gs_data:data[ind_j[m]][ind_i[m]] - mean)) ;
         count++;
     }
     return var/count;
 }



