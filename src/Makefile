ifneq ($(CC),icc)
CXXFLAGS := -g -O3
else
CXX=icpc
CPPFLAGS := -g -fast
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
DEP  := $(SRC:.cpp=.dep)
TDEP := $(TSRC:.cpp=.dep)
BIN  := pet_matric_mc

all: $(BIN)

-include $(DEP) $(TDEP)

$(BIN): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $<

test: $(TOBJ) 
	$(CXX) $(LDFLAGS) -o $@ $<

clean:
	rm -f $(OBJ) $(TOBJ) $(DEP) $(TDEP)
