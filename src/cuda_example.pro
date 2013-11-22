include(common.pri)

CONFIG += cuda
NVCC_FLAGS += -Xptxas -O3 -arch=sm_30

SOURCES += cuda_example/main.cpp
CUDA_SOURCES += cuda_example/kernel.cu
