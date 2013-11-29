TEMPLATE = subdirs

CONFIG += c++11

!CONFIG(verbose): CONFIG += silent

# don't use SDK but call compilers directly on Mac when not using Clang
macx:!macx-clang {
  CONFIG -= sdk
}

SUBDIRS += $$files(src/*.pro)

# tests needs to be explicitely enabled
!test {
  SUBDIRS -= src/test.pro
}

!cuda {
  SUBDIRS -= $$files(src/cuda_*.pro) $$files(src/*_cuda.pro)
}
