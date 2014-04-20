#pragma once

#include <cuda_runtime.h>
#include <cstdio>

/*-----------------------------------------GENERATORS-----------------------------------------------------*/

__constant__ unsigned int shift1[4] = { 6, 2, 13, 3 };
__constant__ unsigned int shift2[4] = { 13, 27, 21, 12 };
__constant__ unsigned int shift3[4] = { 18, 2, 7, 13 };
__constant__ unsigned int offset[4] = { 4294967294,
                                        4294967288,
                                        4294967280,
                                        4294967168 };

__shared__ unsigned int randStates[32];

//-----------end state rng
//---------------------------------------------------------//

__device__ float devData;

__device__ inline unsigned TausStep(unsigned& z,
                                    int S1,
                                    int S2,
                                    int S3,
                                    unsigned M) {
  unsigned b = (((z << S1) ^ z) >> S2);
  return z = (((z & M) << S3) ^ b);
}

__device__ inline unsigned LCGStep(unsigned& z, unsigned A, unsigned C) {
  return z = (A * z + C);
}

__device__ inline float HybridTaus(unsigned& z1,
                                   unsigned& z2,
                                   unsigned& z3,
                                   unsigned& z4) {
  return 2.3283064365387e-10f * (TausStep(z1, 14, 16, 15, 4294967294UL) ^
                                 TausStep(z2, 2, 44, 7, 4294967288UL) ^
                                 TausStep(z3, 3, 144, 17, 4294967280UL) ^
                                 LCGStep(z4, 1664525, 1013904223UL));
}

__device__ float rand_MWC_co(unsigned long long* x, unsigned int* a) {
  // Generate a random number [0,1)
  *x = (*x & 0xffffffffull) * (*a) + (*x >> 32);
  return __fdividef(
      __uint2float_rz((unsigned int)(*x)),
      (float)0x100000000);  // The typecast will truncate the x so that it is
                            // 0<=x<(2^32-1),__uint2float_rz ensures a round
                            // towards zero since 32-bit floating point cannot
                            // represent all integers that large. Dividing by
                            // 2^32 will hence yield [0,1)

}  // end __device__ rand_MWC_co

__device__ float rand_MWC_oc(unsigned long long* x, unsigned int* a) {
  // Generate a random number (0,1]
  return 1.0f - rand_MWC_co(x, a);
}  // end __device__ rand_MWC_oc
