TARGET := $(shell uname -s)

# just in case
ifdef OPENMP
OMP := 1
endif

CXXFLAGS += $(CXXUFLAGS)
LDFLAGS  += $(LDUFLAGS)

ifeq ($(CC),icc)
CXX := icpc
OPT := $(or $(O),fast)
STD := gnu++0x
else
OPT := $(or $(O),3)
ifeq ($(findstring gcc,$(CC)),gcc)
CXX := $(subst gcc,g++,$(CC))
endif
STD := c++11
endif

CXXFLAGS += -std=$(STD)
CPPFLAGS += -g
ifeq ($(OPT),fast)
CPPFLAGS += -fast
else
CPPFLAGS += -O$(OPT)
endif


ifeq ($(TARGET),Darwin)
# force GCC, as clang has not OpenMP on OSX
ifdef OMP
CC  := gcc
CXX := g++
endif
ifneq ($(CXX),g++)
CXXFLAGS += -stdlib=libc++
LDFLAGS  += -stdlib=libc++
endif
endif

ifdef OMP
CPPFLAGS += -fopenmp
LDFLAGS  += -fopenmp
endif

CPPFLAGS += -MMD

# dependencies
CPPFLAGS += -I../lib/catch/include
CPPFLAGS += -I../lib/cmdline

# binary is folder name
BIN  := $(notdir $(realpath .))
SRC  := $(filter-out %_test.cpp, $(wildcard *.cpp))
OBJ  += $(SRC:.cpp=.o)
DEP  += $(SRC:.cpp=.d)

TSRC := $(wildcard *_test.cpp)
TOBJ += t.o $(TSRC:.cpp=.o)
TDEP += t.d $(TSRC:.cpp=.d)

all: $(BIN)

# update submodules
$(OBJ):  ../lib/cmdline/cmdline.h
../lib/cmdline/cmdline.h:
	cd .. && git submodule update --init lib/cmdline
$(TOBJ): ../lib/catch/include/catch.hpp
../lib/catch/include/catch.hpp:
	cd .. && git submodule update --init lib/catch

-include $(DEP) $(TDEP)

# main executable
$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

# unit testing
test: t ; @./t
t: $(TOBJ)
	 $(CXX) $(LDFLAGS) -o $@ $^
t.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ \
	-D CATCH_CONFIG_MAIN \
	-D CATCH_CONFIG_USE_ANSI_COLOUR_CODES \
	-x c++ ../lib/catch/include/catch.hpp

# cleaning
clean:
	rm -f $(OBJ) $(DEP) $(TOBJ) $(TDEP)
