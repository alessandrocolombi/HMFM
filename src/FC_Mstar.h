#ifndef __FC_M_STAR_H__
#define __FC_M_STAR_H__

#include "FullConditional.h"

class FC_Mstar : public FullConditional{

public:
    FC_Mstar(/* args */);
    ~FC_Mstar();
    void update(GS_data& gs_data, sample::GSL_RNG gs_engine) override;
};

#endif