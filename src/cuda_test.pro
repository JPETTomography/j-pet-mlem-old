include(common.pri)

CONFIG += cuda

SOURCES += cuda_test/main.cpp
CUDA_SOURCES += cuda_test/kernel.cu
