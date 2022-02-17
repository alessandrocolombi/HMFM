#ifndef __PARTITION_H__
#define __PARTITION_H__

#include "FullConditional.h"
#include "recurrent_traits.h"

class FullConditional;

class Partition : public FullConditional{
private:
    /* Fix Partition option */
    bool Partition_fixed;
public:
    std::vector< std::vector<unsigned int>> C; // matrix of c_ji (non ordered)
    std::set<unsigned int> s;  // A set which will be useful for ordering components
    std::vector<unsigned int> clust_out; // vector of ORDERED allocated components
    Partition(){};
    Partition(std::string na, const unsigned int d, const std::vector<unsigned int> & n_j, bool FixPart);
    ~Partition(){};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
    double log_norm(double x, double u, double s) const;
};

#endif
