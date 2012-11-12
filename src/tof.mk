VPATH    += ../lib/OpenGL2.0/src/SSG
VPATH    += ../lib/OpenGL2.0/src/OGLUtils
tof_OBJ  += Uniform.o Light.o Shader.o glutils.o glmutils.o logger.o
tof_DEP  += $(tof_OBJ:.o=.d)

CPPFLAGS += -I../lib/OpenGL2.0/src/OGLUtils
CPPFLAGS += -I../lib/OpenGL2.0/src/SSG
CPPFLAGS += -I../lib/glm

ifneq ($(TARGET),Darwin)
tof_LDFLAGS  += -lglut -lGL -lGLEW
tof_LDFLAGS  += -llog4cxx
CPPFLAGS += -DHAVE_LOG4CXX
else
tof_LDFLAGS  += -framework OpenGL
tof_LDFLAGS  += -framework GLUT
endif

# update submodules
$(TOF_OBJ): ../lib/glm/glm/glm.hpp
../lib/glm/glm/glm.hpp:
	cd .. && git submodule update --init lib/glm
$(TOF_OBJ): ../lib/OpenGL2.0/src/OGLUtils/glmutils.h
../lib/OpenGL2.0/src/OGLUtils/glmutils.h:
	cd .. && git submodule update --init lib/OpenGL2.0
