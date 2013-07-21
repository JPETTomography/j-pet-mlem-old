TEMPLATE = subdirs

cache()

equals(PWD, $$OUT_PWD) {
  MAKEFILE = PET.mk
}

SUBDIRS += src/2d_xy_matrix.pro
SUBDIRS += src/2d_xy_phantom.pro
SUBDIRS += src/2d_xy_reconstruction.pro
SUBDIRS += src/2d_strip_reconstruction.pro
SUBDIRS += src/test.pro
