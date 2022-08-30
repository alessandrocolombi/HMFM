#ifndef GDFMM_FC_UMARGINAL_H
#define GDFMM_FC_UMARGINAL_H

#include "FullConditional.h"
#include "FC_U.h"
#include "recurrent_traits.h"

class FC_UMarginal : public FC_U{

public:

    FC_UMarginal(bool _keepfixed):FC_U(_keepfixed){};
    FC_UMarginal(std::string _na, bool _keepfixed):FC_U(_na,_keepfixed){};
    ~FC_UMarginal() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);

};

#endif //GDFMM_FC_U_H
