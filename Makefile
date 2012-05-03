
ifeq "$(strip $(COMPILER))" "gcc"
CC=g++
CXX=g++
CXXFLAGS+=-O3 -g
else
CC=icc -DGTEST_HAS_TR1_TUPLE=0
CXX=icpc -DGTEST_HAS_TR1_TUPLE=0
CXXFLAGS+=-fast -g
endif


LD=$(CXX)
DP=$(CXX) -MM

SED=sed
MV=mv

ifdef OPENMP
ifeq "$(strip $(COMPILER))" "gcc"
CXXFLAGS+=-fopenmp
LDFLAGS +=-fopenmp
else
CXXFLAGS+=-openmp
LDFLAGS +=-openmp
endif
endif


LDFLAGS+=-lboost_program_options

VPATH=$(HOME)/downloads/gtest-1.6.0/src


PREP_C =  perl -pe 's!(${notdir $*})\.o[ :]*!${dir $*}$$1.o $@ : !g' > $@


define make-depend
  $(DP) $(CFLAGS) $(CXXFLAGS) $1 | \
  perl -pe 's!(${notdir $2})[ :]*!${dir $2}$$1 $3 : !g' > $3.tmp
  $(MV) $3.tmp $3
endef

%.o: %.cpp
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) -c $(CXXFLAGS) -o $@ $<






MAIN=pet_matrix_mc.o
#OBJECTS= sweep.o random.o
OBJECTS=
TESTS_OBJECTS=$(subst .cpp,.o,$(wildcard *_test.cpp)) gtest-all.o

DEPEND=$(subst .o,.d,$(MAIN)) $(subst .o,.d,$(OBJECTS)) $(subst .o,.d,$(TESTS_OBJECTS))

-include $(DEPEND)

build: $(MAIN) $(OBJECTS)
	$(LD) $(LDFLAGS) -o pet_matrix_mc  $^	

test: $(OBJECTS) $(TESTS_OBJECTS) 
	$(LD) $(LDFLAGS) -o test  $^  -lpthread

clean:
	rm -f *.o *.d