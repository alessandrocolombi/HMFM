#include "Partition.h"

void Partition::updatePart(GS_data& gs_data, sample::GSL_RNG gs_engine){
  // From gs_data all needed variable are retrived
  unsigned int k = gs_data.K; // number of cluster
  unsigned int d = gs_data.d; // number of group
  std::vector<unsigned int> n_j = gs_data.n_j;//number of observation per group
  unsigned int M = gs_data.M;
  Eigen::Matrix S = gs_data.S;
  std::vector<double> mu=gs_data.mu;
  std::vector<double> sigma=gs_data.sigma;
  std::vector<std::vector<unsigned int>> probs;
  std::set<unsigned int> s;
  std::vector<unsigned int> clust_out;
  double probs_max;

  //Bisogna definire data preso da Gibbs Sampler

  // Genero il vettore probs che mi servirà per il sample
  for(unsigned j=0; j<d; j++){
    std::vector<unsigned int> v(n_j[j]);
    probs.push_back(v); //Create a vector for every J
    for(unsigned i=0; i<n_j[j]; i++){
      for(unsigned m=0; m<M; m++){
        std::normal_distribution<double> d{mu[m],sigma[m]*sigma[m]};
        v.push_back(log(S(j,m) + log(d(data[j][i])); //in every and for every component put the lo likelihood
      }//control if S(j,) is correct
      probs_max=std::max(probs[i])
      for(unsigned m=0; m<M; m++){
        probs[i][m] <- exp(probs[i][m] - probs_max);
        }
    }

    // Assegno tramite il sample su probs a ogni cluster un'etichetta
    //If M==1 populate C matrix with ones
    if (M == 1){
      for (unsigned i=0; i<n_j[j]; i++){
        C[j][i] <- 1;
      }
    }
    else{
      for (unsigned i=0; i<n_j[j]; i++)) {
          std::discrete_distribution<> d(probs[i]); //non so se funziona con i dato che devo accedere
          C[j][i]=std::sample(vec.begin,vec.end(),std::back_inserter(out),1, gs_engine);
          // da mettere anche prob = probs[[j]][i, 1:M],replace = T)
      }
    }/* per ogni dato nel livello j
    Creiamo una matrice della stessa dimensione della matrice dei dati,
    dove ogni riga contiene le etichette non ordinate per ciascun dato di
    quel livello */
  };

  //create vector of allocated components
  for(unsigned j=0; j<d; j++){
    for(unsigned i=0; i<n_j[j]; i++){
      s.insert(C[j][i]);
      clust_out.assign(s.begin(),s.end());
    }
  }
  K=size(clust_out)
  GS_data::set_k(K) //secondo me è necessario un set_K
}
