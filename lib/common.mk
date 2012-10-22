TARGET := $(shell uname -s)

ifneq ($(CC),icc)
CPPFLAGS := -g -O3
else
CXX=icpc
CPPFLAGS := -g -fast
endif

CXXFLAGS += -std=c++11
ifeq ($(TARGET),Darwin)
CXXFLAGS += -stdlib=libc++
LDFLAGS  += -stdlib=libc++
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
$(OBJ): $(SRC)
$(SRC):  ../lib/cmdline/cmdline.h
../lib/cmdline/cmdline.h:
	cd .. && git submodule update --init lib/cmdline
$(TOBJ): $(TSRC)
$(TSRC): ../lib/catch/include/catch.hpp
../lib/catch/include/catch.hpp:
	cd .. && git submodule update --init lib/catch

-include $(DEP) $(TDEP)

# main executable
$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

# unit testing
test: t ; ./t
t: $(TOBJ)
	 $(CXX) $(LDFLAGS) -o $@ $^
t.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ \
	-D CATCH_CONFIG_MAIN -x c++ ../lib/catch/include/catch.hpp

# cleaning
clean:
	rm -f $(OBJ) $(DEP) $(TOBJ) $(TDEP)
