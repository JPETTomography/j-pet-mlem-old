TARGET := $(shell uname -s)

# just in case
ifdef OPENMP
OMP := 1
endif

ifeq ($(CC),icc)
CXX := icpc
OPT := $(or $(O),fast)
else
OPT := $(or $(O),3)
ifeq ($(findstring gcc,$(CC)),gcc)
CXX := $(subst gcc,g++,$(CC))
endif
endif

STD := c++0x

CXXFLAGS += -std=$(STD)
CPPFLAGS += -g
ifeq ($(OPT),fast)
CPPFLAGS += -fast
else
CPPFLAGS += -O$(OPT)
endif

# verbosity (use make Q= to disable)
Q := @

ifeq ($(TARGET),Darwin)
# force GCC, as clang has not OpenMP on OSX
ifdef OMP
ifeq ($(findstring g++,$(CXX))$(findstring icpc,$(CXX)),)
CC  := gcc
CXX := g++
endif
endif
# link to libc++ if not using GCC or ICC
ifeq ($(findstring g++,$(CXX))$(findstring icpc,$(CXX)),)
CXXFLAGS += -stdlib=libc++
LDFLAGS  += -stdlib=libc++
endif
# use static OpenMP
ifneq ($(findstring icpc,$(CXX)),)
LDFLAGS  += -openmp-link=static
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
SRC  := $(filter-out %_test.cpp, $(filter-out pet_%.cpp, $(wildcard *.cpp)))
OBJ  := $(SRC:.cpp=.o)
DEP  += $(SRC:.cpp=.d)

TSRC := $(wildcard *_test.cpp)
TOBJ := t.o $(TSRC:.cpp=.o)
DEP  += t.d $(TSRC:.cpp=.d)

MAIN := $(wildcard pet_*.cpp)
MOBJ := $(MAIN:.cpp=.o)
BIN  := $(patsubst %.cpp,%,$(MAIN))
DEP  += $(MAIN:.cpp=.d)

all: $(BIN)

# update submodules
$(MOBJ) $(OBJ): ../lib/cmdline/cmdline.h
../lib/cmdline/cmdline.h:
	cd .. && git submodule update --init lib/cmdline
$(TOBJ): ../lib/catch/include/catch.hpp
../lib/catch/include/catch.hpp:
	cd .. && git submodule update --init lib/catch

-include $(DEP)
-include $(wildcard *.mk)

# colorful reports
bold   := \033[1m
yellow := \033[33m
cyan   := \033[36m
none   := \033[0m
define CXX_info
	@if [ -t 1 ]; then echo "$(cyan)$(CXX)$(none)\t$<"; else echo "$(CXX) $<"; fi
endef
define LD_info
	@if [ -t 1 ]; then echo "$(yellow)link$(none)\t$(bold)$@$(none)"; else echo "link $@"; fi
endef

%.o: %.cpp
	$(call CXX_info)
	$(Q)$(CXX) $(CXXFLAGS) $(CXXUFLAGS) $(CPPFLAGS) -c -o $@ $<

# main executables
define BIN_template
$(1): $(1).o $$(OBJ) $$($(1)_OBJ)
	$$(call LD_info)
	$(Q)$$(CXX) $$(LDFLAGS) $$(LDUFLAGS) $$($(1)_LDFLAGS) -o $$@ $$^
LOBJ += $$($(1)_OBJ)
endef
$(foreach bin,$(BIN),$(eval $(call BIN_template,$(bin))))

# unit testing
test: t ; @./t
t: $(TOBJ)
	$(call LD_info)
	@$(CXX) $(LDFLAGS) -o $@ $^
t.o:
	@if [ -t 1 ]; then echo "$(green)$(CXX)$(none)\t$<"; else echo "$(CXX) $<"; fi
	$(Q)$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ \
	-D CATCH_CONFIG_MAIN \
	-D CATCH_CONFIG_USE_ANSI_COLOUR_CODES \
	-x c++ ../lib/catch/include/catch.hpp

# cleaning
clean:
	rm -f $(OBJ) $(TOBJ) $(MOBJ) $(LOBJ) $(DEP)
