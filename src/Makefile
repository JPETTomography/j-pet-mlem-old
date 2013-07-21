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
ifneq ($(findstring gcc,$(CC)),)
CXX := $(subst gcc,g++,$(CC))
else
ifneq ($(findstring clang,$(CC)),)
CXX := $(subst clang,clang++,$(CC))
endif
endif
endif

STD := c++11

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

ifdef STATIC
LDFLAGS  += -static
endif

CPPFLAGS += -MMD
CPPFLAGS += -Wall
ifneq ($(findstring g++,$(CXX))$(findstring icpc,$(CXX)),)
# but don't warn about compare sign and comments on GCC
CPPFLAGS += -Wno-sign-compare -Wno-comment
endif


ifeq ($(TARGET),Linux)
ifeq ($(findstring clang,$(CXX)),)
CPPFLAGS += -mssse3 -mfpmath=sse
LDFLAGS += -lrt
endif
endif

# dependencies
CPPFLAGS += -I../lib/catch/include
CPPFLAGS += -I../lib/cmdline
CPPFLAGS += -I.

MAIN := $(wildcard **/*_cmd.cpp)
MOBJ := $(MAIN:.cpp=.o)
BIN  := $(patsubst %_cmd.cpp,%,$(MAIN))
DEP  += $(MAIN:.cpp=.d)

# binary is folder name
SRC  := $(filter-out util/test.cpp, $(filter-out %_test.cpp, $(filter-out %_cmd.cpp, $(wildcard **/*.cpp))))
OBJ  := $(SRC:.cpp=.o)
DEP += $(SRC:.cpp=.d)

TSRC := $(wildcard **/*_test.cpp)
TOBJ := util/test.o $(TSRC:.cpp=.o)
DEP  += util/test.d $(TSRC:.cpp=.d)

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
	@[ -t 1 ] && printf "$(cyan)%s$(none)\t%s\n" "$(CXX)" "$<" || echo "$(CXX) $<"
endef
define LD_info
	@[ -t 1 ] && printf "$(yellow)link$(none)\t$(bold)%s$(none)\n" "$@" || echo "link $@"
endef

%.o: %.cpp
	$(call CXX_info)
	$(Q)$(CXX) $(CXXFLAGS) $(CXXUFLAGS) $(CPPFLAGS) -c -o $@ $<

# main executables
define BIN_template
$(1): $(1)_cmd.o $$(OBJ) $$($(1)_OBJ)
	$$(call LD_info)
	$(Q)$$(CXX) -o $$@ $$^ $$(LDFLAGS) $$(LDUFLAGS) $$($(1)_LDFLAGS)
LOBJ += $$($(1)_OBJ)
endef
$(foreach bin,$(BIN),$(eval $(call BIN_template,$(bin))))

# unit testing
test: util/test ; @util/test
util/test: $(TOBJ) $(OBJ)
	$(call LD_info)
	@$(CXX) -o $@ $^  $(LDFLAGS)

# cleaning
clean:
	rm -f $(OBJ) $(TOBJ) $(MOBJ) $(LOBJ) $(DEP)

style:
	find . \( -iname '*.c' -or -iname '*.cpp' -or -iname '*.h' \) -exec clang-format -style=file -i {} +
