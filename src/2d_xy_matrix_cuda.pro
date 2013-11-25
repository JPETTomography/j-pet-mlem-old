include(common.pri)

CONFIG += cuda
NVCC_FLAGS += -Xptxas -O3 -arch=sm_30

SOURCES += 2d_xy_cuda/matrix_cuda.cpp
CUDA_SOURCES += 2d_xy_cuda/kernel.cu
