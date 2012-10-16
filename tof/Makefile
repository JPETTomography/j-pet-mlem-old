ifneq ($(CC),icc)
CPPFLAGS := -g -O3
else
CXX=icpc
CPPFLAGS := -g -fast
endif

CXXFLAGS += -std=c++11
ifeq ($(shell uname -s),Darwin)
CXXFLAGS += -stdlib=libc++
LDFLAGS  += -stdlib=libc++
endif

ifdef OMP
CPPFLAGS += -fopenmp
LDFLAGS  += -fopenmp
endif

CPPFLAGS += -MMD
LDFLAGS  += -lboost_program_options

CXXFLAGS += -I$(HOME)/downloads/gtest-1.6.0/
CXXFLAGS += -I$(HOME)/PROJECTS/OpenGL2.0/src/SSG/
CXXFLAGS += -I$(HOME)/PROJECTS/OpenGL2.0/src/OGLUtils/

LDFLAGS  += -lglut -llog4cxx -lGLEW

MAIN          := main.o
OBJECTS       := Uniform.o Light.o Shader.o logger.o glutils.o glmutils.o
TESTS_OBJECTS := $($(wildcard *_test.cpp):.cpp=.o) gtest-all.o

DEPEND := $(MAIN:.o=.d) $(OBJECTS:.o=.d) $(TESTS_OBJECTS:.o=.d)

all: tof

-include $(DEPEND)

tof: $(MAIN) $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

test: $(OBJECTS) $(TESTS_OBJECTS)
	$(CXX) $(LDFLAGS) -o test $^ -lpthread

clean:
	rm -f *.o *.d
