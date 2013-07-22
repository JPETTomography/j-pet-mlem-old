TEMPLATE = subdirs

CONFIG += silent

# in Qt5 on Mac we require calling cache of we get complaints
macx {
  greaterThan(QT_MAJOR_VERSION, 4): cache()
}

SUBDIRS += $$files(src/*.pro)
