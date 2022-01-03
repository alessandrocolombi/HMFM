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
        double var(double mean,std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j, const std::vector<std::vector<double>>& data);
        double mean(std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j, const std::vector<std::vector<double>>& data);
    public:
        FC_tau(/* args */);
        ~FC_tau();
        void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
    };


#endif //GDFMM_FC_TAU_H
