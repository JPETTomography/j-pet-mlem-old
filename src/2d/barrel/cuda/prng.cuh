#pragma once

#include <cuda_runtime.h>
#include "util/cuda/compat.h"

__constant__ unsigned int shift1[4] = { 6, 2, 13, 3 };
__constant__ unsigned int shift2[4] = { 13, 27, 21, 12 };
__constant__ unsigned int shift3[4] = { 18, 2, 7, 13 };
__constant__ unsigned int offset[4] = { 4294967294,
                                        4294967288,
                                        4294967280,
                                        4294967168 };

_ inline unsigned TausStep(unsigned& z, int S1, int S2, int S3, unsigned M) {
  unsigned b = (((z << S1) ^ z) >> S2);
  return z = (((z & M) << S3) ^ b);
}

_ inline unsigned LCGStep(unsigned& z, unsigned A, unsigned C) {
  return z = (A * z + C);
}

_ inline float HybridTaus(unsigned& z1,
                          unsigned& z2,
                          unsigned& z3,
                          unsigned& z4) {
  return 2.3283064365387e-10f * (TausStep(z1, 14, 16, 15, 4294967294UL) ^
                                 TausStep(z2, 2, 44, 7, 4294967288UL) ^
                                 TausStep(z3, 3, 144, 17, 4294967280UL) ^
                                 LCGStep(z4, 1664525, 1013904223UL));
}
