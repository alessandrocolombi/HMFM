#include "Partition.h"

void Partition::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
  // From gs_data all needed variable are retrived
  unsigned int k = gs_data.K; // number of cluster
  unsigned int d = gs_data.d; // number of group
  unsigned int M = gs_data.M; // number of components
  const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
  std::vector<unsigned int> n_j = gs_data.n_j;// number of observation per group
  const std::vector<double>& mu = gs_data.mu; // Vector of means
  const std::vector<double>& sigma = gs_data.sigma; // Vector of standard deviations

  std::set<unsigned int> s;// A set which will be useful for cluster
  GDFMM_Traits::VecRow probsvec(M);


  // Define data taken from gs_data
  const std::vector<std::vector<double>>& data = gs_data.data;
  // Initialization of probs_max
  double probs_max;

  //sample::rmultinomial<std::vector<unsigned int>> multinomial;
  sample::sample_index sample_index;
  // std::cout<<d<<std::endl;
  C.clear();
  // Generate matrix of "probabilities" for each observation
  for(unsigned j=0; j<d; j++){
    std::vector<std::vector<double>> probs;
    std::vector<GDFMM_Traits::VecRow> probsmat;
    // matrix of probability we need one of them for every j
    // std::cout<<n_j[j]<<std::endl;
    for(unsigned i=0; i<n_j[j]; i++){
      std::vector<double> vec(M);
      for(unsigned m=0; m<M; m++){
        vec[m]=log(S(j,m)) + log_norm(data[j][i], mu[m], sigma[m]);
        //in every and for every component put the log likelihood
        //std::cout<<data[j][i]-mu[m]<<std::endl;
        //std::cout<<sigma[m]<<std::endl;
        //std::cout<<log(pdfnorm(data[j][i]-mu[m],sigma[m]))<<std::endl;
        //std::cout<<log(S(j,m));
        //std::cout<<log(S(j,m)) + log_norm(data[j][i], mu[m], sigma[m])<<std::endl;
      }

      //std::cout<<std::endl;
      probs.push_back(vec);
      //std::cout<<v[j];
      //Create a vector for eve       //probs è una matrice che ha numero di righe variabile ma sempre M colonne
      probs_max=*max_element(probs[i].begin(), probs[i].end());
      //probs è una matrice che ha numero di righe variabile ma sempre M colonne
      //std::cout<<probs_max<<std::endl;
      for(unsigned m=0; m<M; m++){
        probsvec(m)=exp(probs[i][m] - probs_max);
       // Rcpp::Rcout<<" p:"<<probsvec(m)<<" ";
      //  Rcpp::Rcout<<m;
      }
      //Rcpp::Rcout<<" -- ";
      probsmat.push_back(probsvec);
    }



    // std::cout<<"step 2"<<std::endl;
    // Assegno tramite il sample su probs a ogni cluster un'etichetta
    //If M==1 populate C matrix with ones
    if (M == 1){
      std::vector<unsigned int> v(n_j[j], 1);
      C.push_back(v);
    }
    else{
      std::vector<unsigned int> dis;

      for (unsigned i=0; i<n_j[j]; i++) {
        // VECCHIA VERSIONE
        //double* arrayprobs = &probs[i][0];
        //std::cout << sample_index(gs_engine, probsmat[i])<< "\n";

        dis.push_back(sample_index(gs_engine, probsmat[i]));
        // NUOVA VERSIONE
        // std::vector<unsigned int> sample(M, 0);
        // sample = multinomial(gs_engine, 1, probs[i]);
      }
      C.push_back(dis);
      /* per ogni dato nel livello j
       Creiamo una matrice della stessa dimensione della matrice dei dati,
       dove ogni riga contiene le etichette non ordinate per ciascun dato di
       quel livello */ //
    }


    /*
     for(unsigned i=0; i<n_j[j]; i++){
     for(unsigned m=0; m<M; m++){
     std::cout<<probs[0][m];
     }
     std::cout<<std::endl;
     }
     */
  }
  // std::cout<<"step 3"<<std::endl;
  //create vector of allocated components
  /*for(unsigned int j=0; j<d; j++){
   for(unsigned int i=0; i<n_j[j]; i++){
   std::cout<<C[j][i];
   }
   std::cout<<std::endl;
  }
   */
  // std::cout<<"step 4"<<std::endl;
  clust_out.clear() ; // svuto il vettore clust_out
  for(unsigned int j=0; j<d; j++){
    for(unsigned int i=0; i<n_j[j]; i++){
      s.insert(C[j][i]);
      clust_out.assign(s.begin(),s.end());
    }
  }
  Rcpp::Rcout << "Allocated components:";
  for (auto it = clust_out.begin(); it !=clust_out.end(); ++it)
    Rcpp::Rcout << ' ' << *it;
  Rcpp::Rcout << std::endl;
  k = clust_out.size();

  // std::cout<<"step 5"<<k<<std::endl;
  gs_data.K = k; // updating K in the struct gs_data
  gs_data.initialize_N(k); // initialize N according to new K
  gs_data.update_Ctilde(C, clust_out);
  for(unsigned m=0; m<gs_data.K; m++){
  Rcpp::Rcout<<gs_data.N_k[m]<< " ";
  }
}


double Partition::log_norm(double x, double u, double s) const {
  const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;
  return log((ONE_OVER_SQRT_2PI/std::sqrt(s))*std::exp(-0.5*(x-u)*(x-u)/s));
}
