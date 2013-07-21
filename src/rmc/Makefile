STD := c++0x
CXXFLAGS += -std=$(STD)
CPPFLAGS += -g
CPPFLAGS += -MMD
CPPFLAGS += -Wall

CPPFLAGS += -I../../lib/catch/include
CPPFLAGS += -I../../lib/cmdline

OBJ  := $(SRC:.cpp=.o)
DEP  += $(SRC:.cpp=.d)

TSRC := $(wildcard *_test.cpp)
TOBJ := t.o $(TSRC:.cpp=.o)
DEP  += t.d $(TSRC:.cpp=.d)

$(TOBJ): ../../lib/catch/include/catch.hpp
../../lib/catch/include/catch.hpp:
	cd .. && git submodule update --init lib/catch

-include $(DEP)

# unit testing


test: t ; @./t
t: $(TOBJ) $(OBJ)
	$(call LD_info)
	@$(CXX) -o $@ $^  $(LDFLAGS)

# cleaning
clean:
	rm -f $(OBJ) $(TOBJ) $(MOBJ) $(LOBJ) $(DEP)

style:
	find . \( -iname '*.c' -or -iname '*.cpp' -or -iname '*.h' \) -exec clang-format -style=Chromium -std=c++11 -i {} \;