include(common.pri)

CONFIG += cuda

SOURCES += cuda_example/main.cpp
CUDA_SOURCES += cuda_example/kernel.cu
