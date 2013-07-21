TEMPLATE = subdirs
MAKEFILE = PET.mk

SUBDIRS += 2d_xy/2d_xy_matrix.pro
SUBDIRS += 2d_xy/2d_xy_phantom.pro
SUBDIRS += 2d_xy/2d_xy_reconstruction.pro
SUBDIRS += 2d_strip/2d_strip_reconstruction.pro

HEADERS += geometry/*.h
HEADERS += util/*.h
