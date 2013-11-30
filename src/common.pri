TEMPLATE = app

CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

# fix no symbol in current context on GCC 4.8
QMAKE_CFLAGS_DEBUG += -gdwarf-2

greaterThan(QT_MAJOR_VERSION, 4) {
  CONFIG += object_parallel_to_source
} else:equals(PWD, $$OUT_PWD) {
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

isEmpty(PNGCONFIG) {
  LIBPNGPATHS += /usr/include/png.h \
                 /usr/include/libpng \
                 /usr/local/include/libpng \
                 /opt/X11/include/libpng15 \

  for(path, LIBPNGPATHS):exists($$path) {
    basename = $$basename(path)
    DEFINES += HAVE_LIBPNG
    # only define path if png is not system one
    !equals(basename, png.h) {
      INCLUDEPATH      += $$path
      for(i, 1..2):path = $$dirname(path)
      LIBS             += -L$$path/lib
    }
    LIBS += -lpng
    break()
  }
} else {
  DEFINES        += HAVE_LIBPNG
  QMAKE_CXXFLAGS += $$system($$PNGCONFIG --I_opts)
  LIBS           += $$system($$PNGCONFIG --L_opts --libs) -lz
}

c++11 {
  # workaround for missing old qmake c++11 config
  !greaterThan(QT_MAJOR_VERSION, 4) {
    *-g++-*|*-g++|*-icc|*-icc-*:QMAKE_CXXFLAGS += -std=c++0x
    else:QMAKE_CXXFLAGS += -std=c++11
  } else:*-icc|*-icc-* {
    # c++11 is not enabled in ICC spec
    QMAKE_CXXFLAGS   += -std=c++11
    macx {
      QMAKE_CXXFLAGS += -stdlib=libc++
      QMAKE_LFLAGS   += -stdlib=libc++
    }
  }
}

unix:*-icc|*-icc-* {
  QMAKE_LFLAGS += -static-intel
}
