include(common.pri)

CONFIG  += c++11
CONFIG  += cuda

SOURCES += cuda_example/example.cpp
CUDA_SOURCES += cuda_example/example_kernel.cu
