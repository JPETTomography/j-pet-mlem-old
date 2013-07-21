include(../make/common.pri)

QMAKE_CXXFLAGS += -std=c++11 -fopenmp -msse4.1
QMAKE_LFLAGS   += -fopenmp

SOURCES += reconstruction_cmd.cpp
HEADERS += *.h
