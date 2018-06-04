# Automatically downloads and builds static libpng and zlib libraries
#
# Author: Adam Strzelecki <adam.strzelecki@uj.edu.pl>
#
# Discussion:
#   This script is mainly for Windows & MSVC compatibility which does not ship
#   libpng by default, also there are no well known binary sources.

# overridable settings
set(PNG_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/libpng.ext.dir
  CACHE PATH   "Location on libpng internal installation folder")
set(PNG_URL    "http://downloads.sourceforge.net/project/libpng/libpng15/older-releases/1.5.23/libpng-1.5.23.tar.gz"
  CACHE STRING "URL of libpng source")
set(ZLIB_URL   "http://zlib.net/zlib-1.2.8.tar.gz"
  CACHE STRING "URL of zlib source")

# use function to now pollute global namespace
function(ExtPNG)

  # main dependency
  include(ExternalProject)

  # determine final name for static library
  if(MSVC)
    set(ZLIB_NAME zlibstatic)
    set(PNG_NAME  libpng15_static)
  else()
    set(ZLIB_NAME libz)
    set(PNG_NAME  libpng)
  endif()

  # build zlib externally
  ExternalProject_Add(zlib.ext
    URL ${ZLIB_URL}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/zlib.ext.dir
    INSTALL_DIR ${PNG_PREFIX}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>)

  add_library(zlib UNKNOWN IMPORTED)
  set_property(TARGET zlib PROPERTY IMPORTED_LOCATION
    ${PNG_PREFIX}/lib/${ZLIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
  if(MSVC)
    set_property(TARGET zlib PROPERTY IMPORTED_LOCATION_RELEASE
      ${PNG_PREFIX}/lib/${ZLIB_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
    set_property(TARGET zlib PROPERTY IMPORTED_LOCATION_DEBUG
      ${PNG_PREFIX}/lib/${ZLIB_NAME}d${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
  add_dependencies(zlib zlib.ext)

  # build libpng externally
  ExternalProject_Add(libpng.ext DEPENDS zlib.ext
    URL ${PNG_URL}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/libpng.ext.dir
    INSTALL_DIR ${PNG_PREFIX}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      -DPNG_STATIC:BOOL=ON -DPNG_SHARED:BOOL=OFF)

  add_library(libpng UNKNOWN IMPORTED)
  set_property(TARGET libpng PROPERTY IMPORTED_LOCATION
    ${PNG_PREFIX}/lib/${PNG_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
  if(MSVC)
    set_property(TARGET libpng PROPERTY IMPORTED_LOCATION_RELEASE
      ${PNG_PREFIX}/lib/${PNG_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
    set_property(TARGET libpng PROPERTY IMPORTED_LOCATION_DEBUG
      ${PNG_PREFIX}/lib/${PNG_NAME}d${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
  add_dependencies(libpng zlib libpng.ext)

  # export libraries to parent scope
  set(PNG_LIBRARIES zlib libpng PARENT_SCOPE)
  set(PNG_INCLUDE_DIRS ${PNG_PREFIX}/include PARENT_SCOPE)
  set(PNG_FOUND YES PARENT_SCOPE)
  mark_as_advanced(PNG_LIBRARIES PNG_INCLUDE_DIRS PNG_PREFIX PNG_FOUND)

endfunction()

ExtPNG()
