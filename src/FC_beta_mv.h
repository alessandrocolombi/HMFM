#ifndef GDFMM_FC_BETA_MV_H
#define GDFMM_FC_BETA_MV_H

#include "FullConditional.h"

class FC_beta_mv: public FullConditional {
    private:
        GDFMM_Traits::VecCol beta0; // prior mean
        GDFMM_Traits::MatCol Sigma0; // prior covariance matrix
        GDFMM_Traits::MatCol invSigma0; // prior precision matrix
        GDFMM_Traits::VecCol invSigma0beta0; //Sigma_0^{-1}*beta0

    public:
        FC_beta_mv(std::string na, const GDFMM_Traits::VecCol& _beta0, const GDFMM_Traits::MatCol& _Sigma0, bool _keepfixed) : FullConditional(na,_keepfixed),
                    beta0(_beta0), Sigma0(_Sigma0) {
                        // invert prior covariance matrix
                       GDFMM_Traits::MatCol I( GDFMM_Traits::MatCol::Identity(Sigma0.rows(), Sigma0.cols()) );
                       invSigma0 = Sigma0.llt().solve(I); 
                       invSigma0beta0 = invSigma0 * beta0;
                    };
        ~FC_beta_mv() {};
        void update_data_summary_statistics(GS_data& gs_data);
        void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
};
#endif
