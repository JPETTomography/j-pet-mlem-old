TEMPLATE = app
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += object_parallel_to_source

HEADERS += geometry/*.h
HEADERS += util/*.h

QMAKE_CXXFLAGS += -std=c++11

DESTDIR = ..

INCLUDEPATH += . \
               ../lib/cmdline \
               ../lib/catch/include \

LIBPNGPATHS += /usr/include/libpng12 \
               /usr/include/libpng \
               /usr/local/include/libpng \
               /opt/X11/include/libpng15 \

for(path, LIBPNGPATHS):exists($$path) {
  DEFINES     += HAVE_LIBPNG
  INCLUDEPATH += $$path
  path         = $$dirname(path)
  path         = $$dirname(path)
  LIBS        += -L$$path/lib -lpng
  SOURCES     += util/png_writer.cpp
  break()
}

linux-g++ {
  DEFINES        += OMP=1
  QMAKE_CXXFLAGS += -fopenmp
  LIBS           += -fopenmp
}

macx-clang {
  QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.7
  QMAKE_CXXFLAGS += -stdlib=libc++
  LIBS           += -stdlib=libc++
}
