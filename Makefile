
ifeq "$(strip $(COMPILER))" "gcc"
CC=g++
CXX=g++
CXXFLAGS+=-O3 -g
else
CC=icc -DGTEST_HAS_TR1_TUPLE=0
CXX=icpc -DGTEST_HAS_TR1_TUPLE=0
#CXXFLAGS+=-O1 -g
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
CXXFLAGS+=-I$(HOME)/downloads/gtest-1.6.0/
CXXFLAGS+=-I$(HOME)/PROJECTS/OpenGL2.0/src/SSG/
CXXFLAGS+=-I$(HOME)/PROJECTS/OpenGL2.0/src/OGLUtils/


VPATH=$(HOME)/downloads/gtest-1.6.0/src
VPATH+=$(HOME)/PROJECTS/OpenGL2.0/src/SSG/
VPATH+=$(HOME)/PROJECTS/OpenGL2.0/src/OGLUtils/

PREP_C =  perl -pe 's!(${notdir $*})\.o[ :]*!${dir $*}$$1.o $@ : !g' > $@

LDFLAGS+= -lglut -llog4cxx -lGLEW

define make-depend
  $(DP) $(CFLAGS) $(CXXFLAGS) $1 | \
  perl -pe 's!(${notdir $2})[ :]*!${dir $2}$$1 $3 : !g' > $3.tmp
  $(MV) $3.tmp $3
endef

%.o: %.cpp
	$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) -c $(CXXFLAGS) -o $@ $<






MAIN=main.o
#OBJECTS= sweep.o random.o
OBJECTS=Uniform.o Light.o Shader.o logger.o glutils.o glmutils.o
TESTS_OBJECTS=$(subst .cpp,.o,$(wildcard *_test.cpp)) gtest-all.o

DEPEND=$(subst .o,.d,$(MAIN)) $(subst .o,.d,$(OBJECTS)) $(subst .o,.d,$(TESTS_OBJECTS))

-include $(DEPEND)

tof: $(MAIN) $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@  $^	

test: $(OBJECTS) $(TESTS_OBJECTS) 
	$(LD) $(LDFLAGS) -o test  $^  -lpthread

clean:
	rm -f *.o *.d