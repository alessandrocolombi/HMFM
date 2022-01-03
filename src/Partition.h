#ifndef __PARTITION_H__
#define __PARTITION_H__

#include "FullConditional.h"
#include "recurrent_traits.h"

class FullConditional;

class Partition : public FullConditional{
private:
  /* hyper parameters */

public:
    std::vector<unsigned int> clust_out;
    std::vector< std::vector<unsigned int>> C;
    Partition();
  ~Partition();
  void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
  double normpdf(double x, double u, double s);
};

#endif
