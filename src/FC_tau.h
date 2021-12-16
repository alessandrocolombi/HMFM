//
// Created by ilari on 15/12/2021.
//

#ifndef GDFMM_FC_TAU_H
#define GDFMM_FC_TAU_H

#include "FullConditional.h"
#include "Partition.h"

    class FC_tau: public FullConditional {
    private:
        double nu_0;
        double sigma_0;
        double mu_0;
        double k_0;
        double var(double mean,std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j);
        double mean(std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j);
    public:
        FC_tau(/* args */);
        ~FC_tau();
        void update(GS_data& gs_data, sample::GSL_RNG gs_engine, const string &c, Partition& p) override;



    };


#endif //GDFMM_FC_TAU_H
