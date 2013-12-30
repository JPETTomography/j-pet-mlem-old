CONFIG  += c++11

include(common.pri)

SOURCES += 2d_xy/matrix_cmd.cpp
HEADERS += 2d_xy/*.h \
    2d_xy/TriangleDetector.h

cuda {
  NVCC_FLAGS         += -arch=compute_20 -code=sm_20,sm_30
  NVCC_FLAGS_RELEASE += -Xptxas -O3

  SOURCES      += 2d_xy/cuda/matrix_cuda.cpp
  CUDA_SOURCES += 2d_xy/cuda/matrix_kernel.cu
  HEADERS      += $$files(2d_xy/cuda/*.h) $$files(2d_xy/cuda/*.cuh)
}