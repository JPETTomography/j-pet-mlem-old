#include "simple_lor.h"

int SimpleLOR::n_lors_ = 0;
int SimpleLOR::n_detectors_ = 0;
const SimpleLOR::Comparator SimpleLOR::less;

std::ostream& operator<<(std::ostream& out, const SimpleLOR& lor) {
  return out;
}
