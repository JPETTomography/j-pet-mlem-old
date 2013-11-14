TEMPLATE = app

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

!CONFIG(verbose):CONFIG += silent

# fix no symbol in current context on GCC 4.8
debug:QMAKE_CXXFLAGS += -gdwarf-2
QMAKE_CXXFLAGS += -g -O3 -funroll-loops
LIBS += -g

greaterThan(QT_MAJOR_VERSION, 4) {
  CONFIG += object_parallel_to_source
} else:equals(PWD, $$OUqT_PWD) {
  CONFIG += object_with_source
}

HEADERS += geometry/*.h
HEADERS += util/*.h
HEADERS += math/*.h

SOURCES += util/png_writer.cpp

# drop binaries one level up
equals(PWD, $$OUT_PWD) {
  DESTDIR = ..
} else {
  DESTDIR = $$OUT_PWD/..
}

# don't use SDK but call compilers directly on Mac when not using Clang
macx:!macx-clang {
  CONFIG -= sdk
}

# turn on OpenMP for GCC & ICC
*-g++-*|*-g++|*-icc|*-icc-* {
  QMAKE_CXXFLAGS += -fopenmp
  LIBS           += -fopenmp
}

# enable realtime library for Linux
linux-*:LIBS     += -lrt

INCLUDEPATH += . \
               ../lib/cmdline \
               ../lib/catch/include \

LIBPNGPATHS += /usr/include/png.h \
               /usr/include/libpng \
               /usr/local/include/libpng \
               /opt/X11/include/libpng15 \

for(path, LIBPNGPATHS):exists($$path) {
  basename = $$basename(path)
  DEFINES += HAVE_LIBPNG
  # only define path if png is not system one
  !equals(basename, png.h) {
    INCLUDEPATH += $$path
    path         = $$dirname(path)
    path         = $$dirname(path)
    LIBS        += -L$$path/lib
  }
  LIBS += -lpng
  break()
}

# workaround for missing old qmake c++11 config
CONFIG(c++11):!greaterThan(QT_MAJOR_VERSION, 4) {
  *-g++-*|*-g++|*-icc|*-icc-*:QMAKE_CXXFLAGS += -std=c++0x
  else:QMAKE_CXXFLAGS += -std=c++11
} else {
  # -std=c++11 is not enabled anyway for ICC
  CONFIG(c++11):*-icc|*-icc-*:QMAKE_CXXFLAGS += -std=c++11
}
