ifneq ($(CC),icc)
CPPFLAGS := -g -O3
else
CXX=icpc
CPPFLAGS := -g -fast
endif

PLATFORM := $(shell uname -s)

CXXFLAGS += -std=c++11
ifeq ($(PLATFORM),Darwin)
CXXFLAGS += -stdlib=libc++
LDFLAGS  += -stdlib=libc++
endif

ifdef OMP
CPPFLAGS += -fopenmp
LDFLAGS  += -fopenmp
endif

CPPFLAGS += -MMD
LDFLAGS  += -lboost_program_options

# external libraries
CXXFLAGS += -I../lib/glm
CXXFLAGS += -I../lib/OpenGL2.0/src/SSG
CXXFLAGS += -I../lib/OpenGL2.0/src/OGLUtils
CXXFLAGS += -I../lib/googletest/include
CXXFLAGS += -I../lib/googletest/src
CXXFLAGS += -I../lib/log4cxx/include
LDFLAGS  += -L../lib/log4cxx/src/.libs
LDFLAGS  += -llog4cxx
ifneq ($(PLATFORM),Darwin)
LDFLAGS  += -lglut
else
LDFLAGS  += -lxml2
LDFLAGS  += -framework OpenGL
LDFLAGS  += -framework GLUT
endif

# external sources
VPATH    += ../lib/OpenGL2.0/src/SSG
VPATH    += ../lib/OpenGL2.0/src/OGLUtils
VPATH    += ../lib/googletest/src

MAIN          := main.o
OBJECTS       := Uniform.o Light.o Shader.o glutils.o glmutils.o logger.o
TESTS_OBJECTS := $($(wildcard *_test.cpp):.cpp=.o) gtest-all.o

DEPEND := $(MAIN:.o=.d) $(OBJECTS:.o=.d) $(TESTS_OBJECTS:.o=.d)

all: tof

-include $(DEPEND)

tof: $(MAIN) $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^
ifeq ($(PLATFORM),Darwin)
	install_name_tool -change /System/Library/Frameworks/GLUT.framework/Versions/A/GLUT /Library/Frameworks/GLUT.framework/Versions/A/GLUT $@
endif

test: $(OBJECTS) $(TESTS_OBJECTS)
	$(CXX) $(LDFLAGS) -o test $^ -lpthread

clean:
	rm -f *.o *.d
