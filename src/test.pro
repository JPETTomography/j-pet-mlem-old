include (common.pri)

DIRS = $$files(*)
for(path, DIRS) {
  SOURCES += $$files($$path/*_test.cpp)
}
SOURCES += util/test.cpp
