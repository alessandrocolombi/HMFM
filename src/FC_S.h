#ifndef GDFMM_FC_S_H
#define GDFMM_FC_S_H

#include "FullConditional.h"
#include "recurrent_traits.h"


class FC_S : public FullConditional{
    bool M_fixed;
public:

    FC_S() {name = "S";};
    FC_S(std::string na, bool M_fix): M_fixed(M_fix){name=na;};
    ~FC_S() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

};

#endif //GDFMM_FC_S_H
