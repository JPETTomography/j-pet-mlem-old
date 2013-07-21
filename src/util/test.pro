DIRS = $$files(../*)
for(path, DIRS) {
  SOURCES += $$files($$path/*_test.cpp)
}
