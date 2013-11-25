include(common.pri)

CONFIG += cuda

SOURCES += cuda_example/example.cpp
CUDA_SOURCES += cuda_example/example_kernel.cu
