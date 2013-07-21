TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++11 -fopenmp -msse4.1
QMAKE_LFLAGS += -fopenmp
QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
SOURCES += main.cpp

HEADERS += \
    data_structures.h \
    scintillator.h \
    spet_reconstruction.h \
    bstream.h \
    phantom.h

