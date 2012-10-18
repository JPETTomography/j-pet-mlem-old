VPATH    += ../lib/OpenGL2.0/src/SSG
VPATH    += ../lib/OpenGL2.0/src/OGLUtils
OBJ      += Uniform.o Light.o Shader.o glutils.o glmutils.o logger.o
DEP      += $(OBJ:.o=.d)

-include ../lib/common.mk

CPPFLAGS += -I../lib/OpenGL2.0/src/OGLUtils
CPPFLAGS += -I../lib/OpenGL2.0/src/SSG

ifneq ($(TARGET),Darwin)
LDFLAGS  += -lglut -lGL -lGLEW
else
LDFLAGS  += -lxml2
LDFLAGS  += -framework OpenGL
LDFLAGS  += -framework GLUT
endif

LDFLAGS  += -llog4cxx
