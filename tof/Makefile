-include ../lib/common.mk

ifneq ($(TARGET),Darwin)
LDFLAGS  += -lglut -lGL -lGLEW
else
LDFLAGS  += -lxml2
LDFLAGS  += -framework OpenGL
LDFLAGS  += -framework GLUT
endif
