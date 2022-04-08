#ifndef GDFMM_FC_U_H
#define GDFMM_FC_U_H

#include "FullConditional.h"
#include "recurrent_traits.h"

class FC_U : public FullConditional{

public:

    FC_U(bool _keepfixed){name = "U"; keep_fixed = _keepfixed;};
    FC_U(std::string na, bool _keepfixed){name = na; keep_fixed = _keepfixed;};
    ~FC_U() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;

};

#endif //GDFMM_FC_U_H
