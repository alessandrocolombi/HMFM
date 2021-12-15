#ifndef GDFMM_FC_U_H
#define GDFMM_FC_U_H

#include "FullConditional.h"
#include "recurrent_traits.h"

class FC_U : public FullConditional{

public:
    FC_U(/* args */);
    ~FC_U();
    void update(GS_data& gs_data, sample::GSL_RNG gs_engine) override;

};

#endif //GDFMM_FC_S_H


#endif //GDFMM_FC_U_H
