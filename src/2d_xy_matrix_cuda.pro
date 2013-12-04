include(common.pri)

CONFIG += cuda
NVCC_FLAGS += -arch=compute_20 -code=sm_20,sm_30
NVCC_FLAGS_RELEASE += -Xptxas -O3

SOURCES += 2d_xy_cuda/matrix_cuda.cpp
CUDA_SOURCES += 2d_xy_cuda/kernel.cu

HEADERS += $$files(2d_xy_cuda/*.h) $$files(2d_xy_cuda/*.cuh)
