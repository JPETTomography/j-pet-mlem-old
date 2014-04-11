CONFIG  += c++11

include(common.pri)

SOURCES += 2d_strip/reconstruction_strip_cmd.cpp
HEADERS += 2d_strip/*.h


cuda {
  NVCC_FLAGS         += -arch=compute_30 -code=sm_30
  NVCC_FLAGS_RELEASE += -Xptxas -v -O3 --use_fast_math

  SOURCES      +=
  CUDA_SOURCES += 2d_strip/cuda/reconstruction_strip_gpu.cu
  HEADERS      += $$files(2d_strip/cuda/*.h) $$files(2d_strip/cuda/*.cuh)
}
