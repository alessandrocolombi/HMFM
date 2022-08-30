#ifndef __GDFMM_FC_POSTERIORMARGINAL_HPP__
#define __GDFMM_FC_POSTERIORMARGINAL_HPP__

#include "FullConditional.h"
#include "Partition.h"
#include "recurrent_traits.h"
#include <Rcpp.h>

/*
class FC_PartitionMarginal : public FullConditional{
private:
    // Fix Partition option 
    // Ho un problema. avendo aggiunto in FullConditional il membro keep_fixed, ora questa informazione è ripetuta. Però nel codice viene sempre usata Partition_fixed.
    // Per ora tengo entrambe, sarebbe sa sistemare togliendo Partition_fixed
    bool Partition_fixed;
public:
    std::vector< std::vector<unsigned int>> C; // matrix of c_ji (non ordered)
    std::set<unsigned int> s;  // A set which will be useful for ordering components
    std::vector<unsigned int> clust_out; // vector of ORDERED allocated components
    FC_PartitionMarginal(){};
    FC_PartitionMarginal(std::string na, const unsigned int d, const std::vector<unsigned int> & n_j, bool FixPart);
    ~FC_PartitionMarginal(){};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
    double log_norm(double x, double u, double s) const;
};
*/

class FC_PartitionMarginal : public Partition{
public:	
	FC_PartitionMarginal(){};
    FC_PartitionMarginal(std::string _na, const unsigned int _d, const std::vector<unsigned int>& _n_j, bool _FixPart):Partition(_na,_d,_n_j,_FixPart){};
    ~FC_PartitionMarginal(){};

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);
};

#endif