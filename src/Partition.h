#include "include_headers.h"
#include "recurrent_traits.h"
#include "GS_data.h"
#include "GSL_wrappers.h"
#include "stdio.h"

class Partition{
private:
  /* hyper parameters */

public:
    std::vector<unsigned int> clust_out;
    std::vector< std::vector<unsigned int>> C;
    Partition(std::vector<std::vector<double>> C);
  ~Partition();
  void update(GS_data& gs_data, sample::GSL_RNG gs_engine);
};
