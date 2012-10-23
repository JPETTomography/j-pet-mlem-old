VPATH    += ../lib/OpenGL2.0/src/SSG
VPATH    += ../lib/OpenGL2.0/src/OGLUtils
OBJ      += Uniform.o Light.o Shader.o glutils.o glmutils.o logger.o
DEP      += $(OBJ:.o=.d)

-include ../lib/common.mk

CPPFLAGS += -I../lib/OpenGL2.0/src/OGLUtils
CPPFLAGS += -I../lib/OpenGL2.0/src/SSG
CPPFLAGS += -I../lib/glm

ifneq ($(TARGET),Darwin)
LDFLAGS  += -lglut -lGL -lGLEW
LDFLAGS  += -llog4cxx
CPPFLAGS += -DHAVE_LOG4CXX
else
LDFLAGS  += -framework OpenGL
LDFLAGS  += -framework GLUT
endif

# update submodules
$(OBJ): ../lib/glm/glm/glm.hpp
../lib/glm/glm/glm.hpp:
	cd .. && git submodule update --init lib/glm
$(OBJ): ../lib/OpenGL2.0/src/OGLUtils/glmutils.h
../lib/OpenGL2.0/src/OGLUtils/glmutils.h:
	cd .. && git submodule update --init lib/OpenGL2.0
