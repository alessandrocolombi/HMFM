#ifndef __FC_M_STAR_H__
#define __FC_M_STAR_H__

#include "FullConditional.h"

class FC_Mstar : public FullConditional{
private:
    bool Partition_fixed;
public:

    FC_Mstar(std::string na, bool FixPart, bool _keepfixed) : FullConditional(na,_keepfixed), Partition_fixed(FixPart){};
    FC_Mstar(bool _keepfixed) : FullConditional("Mstar",_keepfixed){};
    ~FC_Mstar() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
};

#endif
