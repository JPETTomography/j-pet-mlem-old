TEMPLATE = subdirs

CONFIG += c++11

!CONFIG(verbose): CONFIG += silent

# don't use SDK but call compilers directly on Mac when not using Clang
macx:!macx-clang {
  CONFIG -= sdk
}

SUBDIRS += $$files(src/*.pro)

# tests & cuda need to be explicitly enabled
!test:SUBDIRS -= src/test.pro
!cuda:SUBDIRS -= src/cuda_example.pro

OTHER_FILES += $$files(scripts/*)
OTHER_FILES += README.md
