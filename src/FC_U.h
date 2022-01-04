#ifndef GDFMM_FC_U_H
#define GDFMM_FC_U_H

#include "FullConditional.h"
#include "recurrent_traits.h"

class FC_U : public FullConditional{

public:
    std::string name="U";
    FC_U(/* args */) {};
    ~FC_U() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

};

#endif //GDFMM_FC_U_H
