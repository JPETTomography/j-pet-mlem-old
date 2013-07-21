TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../../lib/cmdline ../../lib/catch/include

LIBPNGPATHS += /usr/include/libpng
LIBPNGPATHS += /usr/local/include/libpng
LIBPNGPATHS += /opt/X11/include/libpng15
for(path, LIBPNGPATHS):exists($$path) {
  DEFINES     += HAVE_LIBPNG
  INCLUDEPATH += $$path
  path         = $$dirname(path)
  path         = $$dirname(path)
  LIBS        += -L$$path/lib -lpng
  break()
}
