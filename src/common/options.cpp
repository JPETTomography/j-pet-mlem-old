#include "options.h"

namespace Common {

void add_cuda_options(cmdline::parser& cl) {
#if HAVE_CUDA
  cl.add("gpu", 'G', "run on GPU (via CUDA)");
  cl.add<int>("cuda-device", 'D', "CUDA device", cmdline::dontsave, 0);
  cl.add<int>("cuda-blocks", 'B', "CUDA blocks", cmdline::dontsave, 64);
  cl.add<int>(
      "cuda-threads", 'W', "CUDA threads per block", cmdline::dontsave, 512);
#else
  (void)cl;  // mark cl as unsued when not using CUDA
#endif
}

void add_openmp_options(cmdline::parser& cl) {
#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#else
  (void)cl;  // mark cl as unsued when not using OpenMP
#endif
}

}  // Common
