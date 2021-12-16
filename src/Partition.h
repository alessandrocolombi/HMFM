#include "include_headers.h"
#include "recurrent_traits.h"
#include "GS_data.h"
#include "GSL_wrappers.h"
#include "stdio.h"

class Partition{
private:
  /* hyper parameters */
  std::vector< std::vector<double>> C;
public:
  Partition(std::vector<std::vector<double>> C);
  ~Partition();
  void updatePart(GS_data& gs_data, sample::GSL_RNG gs_engine);
};
