TEMPLATE = app

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += object_parallel_to_source
CONFIG += silent

HEADERS += geometry/*.h
HEADERS += util/*.h

QMAKE_CXXFLAGS += -std=c++11

equals(PWD, $$OUT_PWD) {
  DESTDIR = ..
} else {
  DESTDIR = $$OUT_PWD/..
}

INCLUDEPATH += . \
               ../lib/cmdline \
               ../lib/catch/include \

LIBPNGPATHS += /usr/include/png.h \
               /usr/include/libpng \
               /usr/local/include/libpng \
               /opt/X11/include/libpng15 \

for(path, LIBPNGPATHS):exists($$path) {
  basename = $$basename(path)
  DEFINES     += HAVE_LIBPNG
  # only define path if png is not system one
  !equals(basename, png.h) {
    INCLUDEPATH += $$path
    path         = $$dirname(path)
    path         = $$dirname(path)
    LIBS        += -L$$path/lib -lpng
  } else {
    LIBS        += -lpng
  }
  SOURCES     += util/png_writer.cpp
  break()
}

*-g++-* {
  DEFINES        += OMP=1
  QMAKE_CXXFLAGS += -fopenmp
  LIBS           += -fopenmp
}

linux-*:LIBS     += -lrt

macx:QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.7

macx-clang {
  QMAKE_CXXFLAGS += -stdlib=libc++
  LIBS           += -stdlib=libc++
}
