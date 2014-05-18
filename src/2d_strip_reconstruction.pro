CONFIG  += c++11

include(common.pri)

SOURCES += 2d_strip/reconstruction_strip_cmd.cpp
HEADERS += 2d_strip/*.h


cuda {
  NVCC_FLAGS         += -Xptxas -v -arch=compute_30 -code=sm_30 -ftz=true -prec-div=false -prec-sqrt=false
  NVCC_FLAGS_RELEASE +=  -O3

  SOURCES      +=
  CUDA_SOURCES += 2d_strip/cuda/reconstruction_strip_gpu.cu
  HEADERS      += $$files(2d_strip/cuda/*.h) $$files(2d_strip/cuda/*.cuh)
}
