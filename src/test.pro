include (common.pri)

CONFIG  += c++11

DIRS = $$files(*)
for(path, DIRS) {
  SOURCES += $$files($$path/*_test.cpp)
}
SOURCES += util/test.cpp
