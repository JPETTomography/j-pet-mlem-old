CONFIG  += c++11
CONFIG  += cuda

include(common.pri)

SOURCES += cuda_example/example.cpp
CUDA_SOURCES += cuda_example/example_kernel.cu
