TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../lib/cmdline ../lib/catch/include

mac:QMAKE_CXXFLAGS += -stdlib=libc++ -std=c++0x
mac:QMAKE_LFLAGS   += -stdlib=libc++
mac:QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.7
mac:INCLUDEPATH += /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/c++/v1
mac:INCLUDEPATH += /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/include
