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

SRC  := $(wildcard *.cpp)
TSRC := $(filter %_test.cpp, $(SRC))
SRC  := $(filter-out %_test.cpp, $(SRC))
OBJ  := $(SRC:.cpp=.o)
TOBJ := $(TSRC:.cpp=.o)
DEP  := $(SRC:.cpp=.d)
TDEP := $(TSRC:.cpp=.d)
BIN  := pet_matrix_mc

all: $(BIN)

-include $(DEP) $(TDEP)

$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $<

test: $(TOBJ) 
	$(CXX) $(LDFLAGS) -o $@ $<

clean:
	rm -f $(OBJ) $(TOBJ) $(DEP) $(TDEP)
