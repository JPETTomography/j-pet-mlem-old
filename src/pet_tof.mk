VPATH    += ../lib/OpenGL2.0/src/SSG
VPATH    += ../lib/OpenGL2.0/src/OGLUtils

pet_tof_OBJ += Uniform.o Light.o Shader.o glutils.o glmutils.o logger.o
pet_tof_DEP += $(tof_OBJ:.o=.d)

CPPFLAGS += -I../lib/OpenGL2.0/src/OGLUtils
CPPFLAGS += -I../lib/OpenGL2.0/src/SSG
CPPFLAGS += -I../lib/glm

ifneq ($(TARGET),Darwin)
pet_tof_LDFLAGS  += -lglut -lGL -lGLEW
pet_tof_LDFLAGS  += -llog4cxx
CPPFLAGS += -DHAVE_LOG4CXX
else
pet_tof_LDFLAGS  += -framework OpenGL
pet_tof_LDFLAGS  += -framework GLUT
endif

# update submodules
pet_tof.o: ../lib/glm/glm/glm.hpp
../lib/glm/glm/glm.hpp:
	cd .. && git submodule update --init lib/glm
pet_tof.o: ../lib/OpenGL2.0/src/OGLUtils/glmutils.h
../lib/OpenGL2.0/src/OGLUtils/glmutils.h:
	cd .. && git submodule update --init lib/OpenGL2.0
